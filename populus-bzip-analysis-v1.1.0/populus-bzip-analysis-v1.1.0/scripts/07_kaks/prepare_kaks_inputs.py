#!/usr/bin/env python3
"""Prepare clean Ka/Ks inputs for the revised 79-OGG Populus bZIP set."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    record_id: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if record_id is not None:
                    if record_id in records:
                        raise ValueError(f"Duplicate FASTA ID in {path}: {record_id}")
                    records[record_id] = "".join(sequence).upper()
                record_id = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)
    if record_id is not None:
        if record_id in records:
            raise ValueError(f"Duplicate FASTA ID in {path}: {record_id}")
        records[record_id] = "".join(sequence).upper()
    return records


def write_fasta(records: dict[str, str], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record_id, sequence in records.items():
            handle.write(f">{record_id}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index:index + 60] + "\n")


def translate(cds: str) -> str:
    return "".join(CODON_TABLE.get(cds[index:index + 3], "X") for index in range(0, len(cds) - 2, 3))


def trim_terminal_stop(cds: str) -> str:
    cds = cds.upper().replace("U", "T")
    return cds[:-3] if translate(cds).endswith("*") else cds


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_orthogroups(path: Path) -> tuple[list[str], dict[str, list[tuple[str, str]]]]:
    groups: dict[str, list[tuple[str, str]]] = {}
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        species_columns = [column for column in (reader.fieldnames or []) if column != "Orthogroup"]
        for row in reader:
            members: list[tuple[str, str]] = []
            for species in species_columns:
                value = row.get(species, "")
                if value:
                    members.extend((species, item.strip()) for item in value.split(",") if item.strip())
            groups[row["Orthogroup"]] = members
    return species_columns, groups


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--nomenclature", type=Path, required=True)
    parser.add_argument("--subfamily-map", type=Path, required=True)
    parser.add_argument("--current-protein-dir", type=Path, required=True)
    parser.add_argument("--legacy-cds", type=Path, required=True)
    parser.add_argument("--retained-cds", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()

    if args.output_root.exists():
        raise SystemExit(f"Refusing to overwrite existing output root: {args.output_root}")
    input_dir = args.output_root / "input"
    pair_dir = args.output_root / "pairs"
    log_dir = args.output_root / "logs"
    source_dir = args.output_root / "source_snapshot"
    for directory in (input_dir, pair_dir, log_dir, source_dir):
        directory.mkdir(parents=True, exist_ok=False)

    species_columns, groups = parse_orthogroups(args.orthogroups)
    nomenclature_rows = read_table(args.nomenclature)
    og_to_nom = {row["New_Orthogroup"]: row for row in nomenclature_rows}
    subfamily_rows = read_table(args.subfamily_map)
    ptbzip_to_subfamily = {row["PtbZIP_ID"]: row["subfamily"] for row in subfamily_rows}

    current_proteins: dict[str, str] = {}
    for path in sorted(args.current_protein_dir.glob("*.fa")):
        for record_id, sequence in read_fasta(path).items():
            if record_id in current_proteins:
                raise ValueError(f"Duplicate current protein ID: {record_id}")
            current_proteins[record_id] = sequence
    legacy_cds = read_fasta(args.legacy_cds)
    retained_cds = read_fasta(args.retained_cds)

    member_to_og: dict[str, str] = {}
    member_to_species: dict[str, str] = {}
    for orthogroup, members in groups.items():
        if orthogroup not in og_to_nom:
            raise ValueError(f"Missing nomenclature for {orthogroup}")
        for species, record_id in members:
            if record_id in member_to_og:
                raise ValueError(f"Member assigned more than once: {record_id}")
            member_to_og[record_id] = orthogroup
            member_to_species[record_id] = species

    revised_cds: dict[str, str] = {}
    revised_peptide: dict[str, str] = {}
    qc_rows: list[dict[str, object]] = []
    for record_id in sorted(member_to_og):
        current_peptide = current_proteins.get(record_id, "").rstrip("*")
        if not current_peptide:
            raise ValueError(f"Current protein missing for assigned member: {record_id}")
        if record_id in retained_cds:
            cds_source = "audited_retained_model"
            cds = trim_terminal_stop(retained_cds[record_id])
        elif record_id in legacy_cds:
            cds_source = "annotated_model_cds"
            cds = trim_terminal_stop(legacy_cds[record_id])
        else:
            raise ValueError(f"CDS missing for assigned member: {record_id}")

        source_translation = translate(cds)
        if "*" in source_translation or "X" in source_translation:
            raise ValueError(f"Invalid CDS translation for {record_id}")

        if current_peptide == source_translation:
            compatibility = "exact"
            final_cds = cds
            final_peptide = current_peptide
            alignment_peptide_origin = "current_orthofinder_protein"
        else:
            offset = source_translation.find(current_peptide)
            if offset >= 0:
                compatibility = "current_is_exact_subsequence_of_cds_translation"
                final_cds = cds[offset * 3:(offset + len(current_peptide)) * 3]
                final_peptide = current_peptide
                alignment_peptide_origin = "current_orthofinder_protein_codon_cropped"
            elif current_peptide.find(source_translation) >= 0:
                compatibility = "cds_translation_is_exact_subsequence_of_current"
                final_cds = cds
                final_peptide = source_translation
                alignment_peptide_origin = "cds_translation_shorter_than_current_protein"
            else:
                compatibility = "not_exactly_compatible"
                final_cds = cds
                final_peptide = source_translation
                alignment_peptide_origin = "cds_translation_for_incompatible_annotation_pair"

        if len(final_cds) % 3 or translate(final_cds) != final_peptide:
            raise ValueError(f"Final CDS-peptide mismatch for {record_id}")
        revised_cds[record_id] = final_cds
        revised_peptide[record_id] = final_peptide

        orthogroup = member_to_og[record_id]
        nom = og_to_nom[orthogroup]
        pangene = nom["Pangene_ID"]
        qc_rows.append({
            "Record_ID": record_id,
            "Species": member_to_species[record_id],
            "Orthogroup": orthogroup,
            "Pangene_ID": pangene,
            "Pangenome_class": nom["Category"],
            "Subfamily": ptbzip_to_subfamily.get(pangene, "NA"),
            "Source_status": "supplementary_model" if record_id in retained_cds else "annotation_derived",
            "CDS_source": cds_source,
            "Current_protein_length_aa": len(current_peptide),
            "KaKs_alignment_peptide_length_aa": len(final_peptide),
            "KaKs_CDS_length_bp": len(final_cds),
            "CDS_protein_compatibility": compatibility,
            "Alignment_peptide_origin": alignment_peptide_origin,
        })

    write_fasta(revised_cds, input_dir / "Populus_bZIP_1762_assigned.current79.cds.fasta")
    write_fasta(revised_peptide, input_dir / "Populus_bZIP_1762_assigned.current79.protein.fasta")

    with (input_dir / "Populus_bZIP_1762_assigned.current79.sequence_qc.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(qc_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(qc_rows)

    pair_rows: list[dict[str, object]] = []
    seen_unordered: set[tuple[str, str]] = set()
    anchor_rows: list[dict[str, object]] = []
    for orthogroup, members in groups.items():
        if len(members) <= 1:
            continue
        ptr_anchors = [record_id for species, record_id in members if species == "Populus_trichocarpa"]
        if ptr_anchors:
            anchors = sorted(ptr_anchors)
            anchor_rule = "all_Populus_trichocarpa_members"
        else:
            max_length = max(len(revised_peptide[record_id]) for _, record_id in members)
            anchors = [min(record_id for _, record_id in members if len(revised_peptide[record_id]) == max_length)]
            anchor_rule = "longest_CDS_compatible_protein_no_Populus_trichocarpa"
        nom = og_to_nom[orthogroup]
        pangene = nom["Pangene_ID"]
        for anchor in anchors:
            anchor_rows.append({
                "Orthogroup": orthogroup,
                "Pangene_ID": pangene,
                "Anchor_gene": anchor,
                "Anchor_species": member_to_species[anchor],
                "Anchor_protein_length_aa": len(revised_peptide[anchor]),
                "Anchor_rule": anchor_rule,
            })
            for _, partner in members:
                if partner == anchor:
                    continue
                unordered = tuple(sorted((anchor, partner)))
                if unordered in seen_unordered:
                    continue
                seen_unordered.add(unordered)
                pair_rows.append({
                    "Pair_ID": f"pair{len(pair_rows) + 1:06d}",
                    "Orthogroup": orthogroup,
                    "Pangene_ID": pangene,
                    "Pangenome_class": nom["Category"],
                    "Subfamily": ptbzip_to_subfamily.get(pangene, "NA"),
                    "Gene1": anchor,
                    "Gene2": partner,
                    "Gene1_species": member_to_species[anchor],
                    "Gene2_species": member_to_species[partner],
                    "Gene1_source_status": "supplementary_model" if anchor in retained_cds else "annotation_derived",
                    "Gene2_source_status": "supplementary_model" if partner in retained_cds else "annotation_derived",
                    "Anchor_rule": anchor_rule,
                })

    with (pair_dir / "current79_homolog_pairs.2col.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        for row in pair_rows:
            handle.write(f"{row['Gene1']}\t{row['Gene2']}\n")
    with (pair_dir / "current79_homolog_pairs.annotated.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(pair_rows)
    with (pair_dir / "current79_anchor_selection.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(anchor_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(anchor_rows)

    source_paths = [args.orthogroups, args.nomenclature, args.subfamily_map, args.legacy_cds, args.retained_cds]
    with (source_dir / "source_files.sha256.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Path", "SHA256"])
        for path in source_paths:
            writer.writerow([str(path), sha256(path)])

    summary = {
        "species": len(species_columns),
        "orthogroups": len(groups),
        "assigned_members": len(member_to_og),
        "input_cds": len(revised_cds),
        "input_proteins": len(revised_peptide),
        "pairs": len(pair_rows),
        "anchors": len(anchor_rows),
        "cds_protein_compatibility": dict(Counter(str(row["CDS_protein_compatibility"]) for row in qc_rows)),
        "alignment_peptide_origins": dict(Counter(str(row["Alignment_peptide_origin"]) for row in qc_rows)),
        "pair_anchor_rules": dict(Counter(str(row["Anchor_rule"]) for row in pair_rows)),
        "unordered_pair_deduplication": True,
    }
    (args.output_root / "PREPARATION_SUMMARY.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
