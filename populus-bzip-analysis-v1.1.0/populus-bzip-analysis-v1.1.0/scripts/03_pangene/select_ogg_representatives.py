#!/usr/bin/env python3
"""Select revised OGG representatives and compare them with historical representatives."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--old-representatives-tsv", required=True, type=Path)
    parser.add_argument("--old-representatives-fasta", required=True, type=Path)
    parser.add_argument("--old-to-new-mapping", required=True, type=Path)
    parser.add_argument("--unassigned", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    if identifier in records:
                        raise RuntimeError(f"Duplicate FASTA ID {identifier} in {path}")
                    records[identifier] = "".join(chunks).upper()
                identifier = line[1:].split()[0]
                chunks = []
            else:
                if identifier is None:
                    raise RuntimeError(f"Sequence before header in {path}")
                chunks.append(line)
    if identifier is not None:
        if identifier in records:
            raise RuntimeError(f"Duplicate FASTA ID {identifier} in {path}")
        records[identifier] = "".join(chunks).upper()
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def split_semicolon(value: str) -> list[str]:
    return [item for item in value.split(";") if item and item != "NA"]


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite output: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    ogg_rows = read_tsv(args.orthogroups)
    if not ogg_rows:
        raise RuntimeError("Orthogroups table is empty")
    species_columns = [field for field in ogg_rows[0] if field != "Orthogroup"]
    if len(species_columns) != 19:
        raise RuntimeError(f"Expected 19 species columns, found {len(species_columns)}")

    all_sequences: dict[str, str] = {}
    sequence_file_by_id: dict[str, str] = {}
    for fasta in sorted(args.input_dir.glob("*.fa")):
        for identifier, sequence in read_fasta(fasta).items():
            if identifier in all_sequences:
                raise RuntimeError(f"Duplicate protein ID across input FASTAs: {identifier}")
            all_sequences[identifier] = sequence
            sequence_file_by_id[identifier] = fasta.name
    if len(all_sequences) != 1764:
        raise RuntimeError(f"Expected 1764 revised input proteins, found {len(all_sequences)}")

    class_rows = read_tsv(args.classification)
    class_by_ogg = {row["Orthogroup"]: row for row in class_rows}
    old_rows = read_tsv(args.old_representatives_tsv)
    old_by_ogg = {row["Orthogroup"]: row for row in old_rows}
    old_fasta_raw = read_fasta(args.old_representatives_fasta)
    old_sequence_by_ogg: dict[str, str] = {}
    old_member_by_ogg: dict[str, str] = {}
    for header, sequence in old_fasta_raw.items():
        old_ogg, separator, member = header.partition("__")
        if not separator:
            raise RuntimeError(f"Unexpected old representative header: {header}")
        old_sequence_by_ogg[old_ogg] = sequence
        old_member_by_ogg[old_ogg] = member
    if set(old_by_ogg) != set(old_sequence_by_ogg):
        raise RuntimeError("Historical representative TSV and FASTA OGG sets differ")

    new_rows: list[dict[str, object]] = []
    representative_fasta: list[tuple[str, str]] = []
    new_rep_member: dict[str, str] = {}
    new_rep_sequence: dict[str, str] = {}
    members_seen: set[str] = set()
    max_length_ties = 0

    for ogg_row in ogg_rows:
        ogg = ogg_row["Orthogroup"]
        if ogg not in class_by_ogg:
            raise RuntimeError(f"Missing classification row for {ogg}")
        class_row = class_by_ogg[ogg]
        members_by_species: dict[str, list[str]] = {}
        for species in species_columns:
            cell = ogg_row.get(species, "").strip()
            members = [item.strip() for item in cell.split(",") if item.strip()]
            members_by_species[species] = members
            for member in members:
                if member not in all_sequences:
                    raise RuntimeError(f"{ogg} member absent from revised FASTAs: {member}")
                if member in members_seen:
                    raise RuntimeError(f"Protein assigned to more than one OGG: {member}")
                members_seen.add(member)

        all_members = [member for species in species_columns for member in members_by_species[species]]
        ptr_members = members_by_species["Populus_trichocarpa"]
        pool = ptr_members if ptr_members else all_members
        rule = "longest_Ptr_in_OGG" if ptr_members else "longest_in_OGG_no_Ptr"
        max_length = max(len(all_sequences[member]) for member in pool)
        tied = sorted(member for member in pool if len(all_sequences[member]) == max_length)
        selected = tied[0]
        selected_species = next(
            species for species, members in members_by_species.items() if selected in members
        )
        selected_gene = selected.rsplit("__", 1)[-1]
        selected_sequence = all_sequences[selected]
        representative_header = f"{ogg}__{selected}"
        contributors = split_semicolon(class_row.get("Old_orthogroup_contributors", ""))
        contributor_pangenes = split_semicolon(class_row.get("Old_pangene_contributors", ""))
        id_matches = [
            old_ogg
            for old_ogg in contributors
            if old_by_ogg.get(old_ogg, {}).get("Representative_ID") == selected
        ]
        sequence_matches = [
            old_ogg
            for old_ogg in contributors
            if old_sequence_by_ogg.get(old_ogg) == selected_sequence
        ]
        if not contributors:
            comparison = "new_OGG_without_old_contributor"
        elif len(contributors) == 1:
            if id_matches:
                comparison = "unchanged_ID_and_sequence"
            elif sequence_matches:
                comparison = "same_sequence_different_ID"
            else:
                comparison = "changed_sequence"
        elif sequence_matches:
            comparison = "matches_one_or_more_merged_old_representatives"
        else:
            comparison = "changed_sequence_after_merge"

        if len(tied) > 1:
            max_length_ties += 1
        source_type = "audited_model" if selected_gene.endswith(".final") else "annotated"
        new_row: dict[str, object] = {
            "New_Orthogroup": ogg,
            "Category": class_row["Category"],
            "Species_present": class_row["Species_present"],
            "Total_genes": class_row["Total_genes"],
            "Representative_member_ID": selected,
            "Representative_header": representative_header,
            "Species": selected_species,
            "Original_gene": selected_gene,
            "Source_type": source_type,
            "Sequence_file": sequence_file_by_id[selected],
            "Seq_length": len(selected_sequence),
            "Selection_rule": rule,
            "Max_length_tie_count": len(tied),
            "Tie_break_rule": "lexicographically_first_ID" if len(tied) > 1 else "NA",
            "Old_orthogroup_contributors": ";".join(contributors) or "NA",
            "Old_pangene_contributors": ";".join(contributor_pangenes) or "NA",
            "Merge_of_multiple_old_OGGs": class_row.get("Merge_of_multiple_old_OGGs", "NA"),
            "Matching_old_representative_ID_OGGs": ";".join(id_matches) or "NA",
            "Matching_old_representative_sequence_OGGs": ";".join(sequence_matches) or "NA",
            "Representative_comparison": comparison,
        }
        new_rows.append(new_row)
        representative_fasta.append((representative_header, selected_sequence))
        new_rep_member[ogg] = selected
        new_rep_sequence[ogg] = selected_sequence

    if len(new_rows) != 79:
        raise RuntimeError(f"Expected 79 revised OGG representatives, found {len(new_rows)}")
    if len(members_seen) != 1762:
        raise RuntimeError(f"Expected 1762 assigned proteins, found {len(members_seen)}")

    mapping_rows = [
        row
        for row in read_tsv(args.old_to_new_mapping)
        if row.get("dataset") == "revised_audited_25"
    ]
    mapping_by_old = {row["Old_Orthogroup"]: row for row in mapping_rows}
    old_comparison_rows: list[dict[str, object]] = []
    for old_ogg in sorted(old_by_ogg):
        old_row = old_by_ogg[old_ogg]
        mapping = mapping_by_old.get(old_ogg)
        mapped_new = (
            split_semicolon(mapping.get("Mapped_new_OGGs", "")) if mapping is not None else []
        )
        mapped_new_reps = [new_rep_member[new_ogg] for new_ogg in mapped_new]
        old_member = old_by_ogg[old_ogg]["Representative_ID"]
        old_sequence = old_sequence_by_ogg[old_ogg]
        id_retained = [new_ogg for new_ogg in mapped_new if new_rep_member[new_ogg] == old_member]
        sequence_retained = [
            new_ogg for new_ogg in mapped_new if new_rep_sequence[new_ogg] == old_sequence
        ]
        mapped_new_is_merge = any(
            len(split_semicolon(class_by_ogg[new_ogg].get("Old_orthogroup_contributors", ""))) > 1
            for new_ogg in mapped_new
        )
        if not mapped_new:
            status = "old_OGG_lost"
        elif len(mapped_new) > 1:
            status = (
                "old_representative_retained_in_split"
                if sequence_retained
                else "old_representative_not_retained_after_split"
            )
        elif mapped_new_is_merge:
            status = (
                "old_representative_retained_in_merged_new_OGG"
                if sequence_retained
                else "old_representative_not_retained_after_merge"
            )
        elif id_retained:
            status = "unchanged_ID_and_sequence_one_to_one"
        elif sequence_retained:
            status = "same_sequence_different_ID_one_to_one"
        else:
            status = "changed_sequence_one_to_one"

        old_comparison_rows.append(
            {
                "Old_Orthogroup": old_ogg,
                "Old_pangene_ID": old_row.get("Pangene_ID", "NA"),
                "Old_representative_ID": old_member,
                "Old_representative_length": len(old_sequence),
                "Mapped_new_OGGs": ";".join(mapped_new) or "NA",
                "Mapped_new_representative_IDs": ";".join(mapped_new_reps) or "NA",
                "Old_representative_ID_retained_in": ";".join(id_retained) or "NA",
                "Old_representative_sequence_retained_in": ";".join(sequence_retained) or "NA",
                "Mapping_status": mapping.get("Mapping_status", "unmapped") if mapping else "unmapped",
                "Representative_status": status,
            }
        )

    unassigned_rows = read_tsv(args.unassigned)
    unassigned_ids: list[str] = []
    for row in unassigned_rows:
        for species in species_columns:
            unassigned_ids.extend(
                item.strip() for item in row.get(species, "").split(",") if item.strip()
            )
    if len(unassigned_ids) != 2:
        raise RuntimeError(f"Expected 2 unassigned proteins, found {len(unassigned_ids)}")

    new_fields = list(new_rows[0])
    old_fields = list(old_comparison_rows[0])
    write_tsv(args.output_dir / "New79_OGG_representatives.tsv", new_rows, new_fields)
    write_fasta(
        args.output_dir / "New79_OGG_representatives.full_length.fasta",
        representative_fasta,
    )
    write_tsv(
        args.output_dir / "Old86_to_New79_representative_comparison.tsv",
        old_comparison_rows,
        old_fields,
    )
    (args.output_dir / "Unassigned_proteins_excluded_from_representative_tree.txt").write_text(
        "\n".join(unassigned_ids) + "\n", encoding="utf-8"
    )

    new_status_counts = Counter(str(row["Representative_comparison"]) for row in new_rows)
    old_status_counts = Counter(str(row["Representative_status"]) for row in old_comparison_rows)
    summary = {
        "selection_rule": {
            "primary": "longest P. trichocarpa protein in each OGG",
            "fallback": "longest protein in OGG when P. trichocarpa is absent",
            "tie_break": "lexicographically first full member ID",
        },
        "revised_input_proteins": len(all_sequences),
        "assigned_proteins": len(members_seen),
        "unassigned_proteins_excluded": len(unassigned_ids),
        "revised_OGGs": len(new_rows),
        "representatives_selected": len(representative_fasta),
        "representatives_from_P_trichocarpa": sum(
            row["Species"] == "Populus_trichocarpa" for row in new_rows
        ),
        "representatives_from_other_species": sum(
            row["Species"] != "Populus_trichocarpa" for row in new_rows
        ),
        "representatives_from_audited_models": sum(
            row["Source_type"] == "audited_model" for row in new_rows
        ),
        "OGGs_with_max_length_ties": max_length_ties,
        "new_OGG_representative_comparison_counts": dict(sorted(new_status_counts.items())),
        "old_OGG_representative_status_counts": dict(sorted(old_status_counts.items())),
        "all_internal_checks_passed": True,
    }
    (args.output_dir / "representative_selection_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
