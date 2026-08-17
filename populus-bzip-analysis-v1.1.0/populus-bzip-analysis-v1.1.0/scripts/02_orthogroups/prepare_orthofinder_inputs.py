#!/usr/bin/env python3
"""Build audited and sensitivity OrthoFinder inputs for the Populus bZIP revision.

The original 52 BITACORA-derived candidates are removed from every data set.
The audited data set receives all 25 validated final models, the
no-independent-new sensitivity set receives only the four existing-locus
annotation corrections, and the strict annotated-only sensitivity set receives
no candidate-derived model.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
from collections import Counter, defaultdict
from pathlib import Path


DATASETS = {
    "revised_audited_25": "all_retained",
    "sensitivity_no_independent_new": "annotation_corrections_only",
    "sensitivity_annotated_only": "annotated_only",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--old-root", required=True, type=Path)
    parser.add_argument("--audit-tsv", required=True, type=Path)
    parser.add_argument("--retained-proteins", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> list[tuple[str, str, str]]:
    records: list[tuple[str, str, str]] = []
    header: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header.split()[0], header, "".join(parts)))
                header = line[1:].strip()
                parts = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                parts.append(line)
    if header is not None:
        records.append((header.split()[0], header, "".join(parts)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")


def safe_species(species: str) -> str:
    return re.sub(r"[^A-Za-z0-9_]", "_", species)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    old_input = args.old_root / "OrthoFinder_input"
    old_mapping_path = args.old_root / "tables" / "bZIP_orthofinder_id_mapping.tsv"

    for path in (old_input, old_mapping_path, args.audit_tsv, args.retained_proteins):
        if not path.exists():
            raise FileNotFoundError(path)
    if args.output_root.exists():
        raise FileExistsError(
            f"Refusing to overwrite existing revision directory: {args.output_root}"
        )

    audit_rows = read_tsv(args.audit_tsv)
    if len(audit_rows) != 52:
        raise RuntimeError(f"Expected 52 audit records, found {len(audit_rows)}")
    retained_rows = [
        row for row in audit_rows if row["retained_in_revised_bZIP_set"] == "Yes"
    ]
    corrections = [
        row for row in retained_rows if row["independent_new_gene"] == "No"
    ]
    independent = [
        row for row in retained_rows if row["independent_new_gene"] == "Yes"
    ]
    if (len(retained_rows), len(independent), len(corrections)) != (25, 21, 4):
        raise RuntimeError(
            "Audit totals do not match 25 retained = 21 independent + 4 corrections: "
            f"{len(retained_rows)}, {len(independent)}, {len(corrections)}"
        )

    old_mapping = read_tsv(old_mapping_path)
    old_new_rows = [row for row in old_mapping if row["Type"] == "new_predicted"]
    old_annotated_rows = [row for row in old_mapping if row["Type"] == "annotated"]
    if (len(old_new_rows), len(old_annotated_rows)) != (52, 1739):
        raise RuntimeError(
            f"Old mapping expected 52 new + 1739 annotated, observed "
            f"{len(old_new_rows)} + {len(old_annotated_rows)}"
        )
    old_new_ids = {row["OrthoFinder_ID"] for row in old_new_rows}
    audit_keys = {(row["species"], row["candidate_gene_id"]) for row in audit_rows}
    mapping_keys = {(row["Species"], row["Original_ID"]) for row in old_new_rows}
    if audit_keys != mapping_keys:
        missing_audit = sorted(mapping_keys - audit_keys)
        missing_mapping = sorted(audit_keys - mapping_keys)
        raise RuntimeError(
            f"Audit/mapping candidate mismatch; missing audit={missing_audit}; "
            f"missing mapping={missing_mapping}"
        )

    old_records_by_species: dict[str, list[tuple[str, str]]] = {}
    all_old_ids: list[str] = []
    for fasta in sorted(old_input.glob("*.fa")):
        species_file = fasta.stem
        records = [(identifier, sequence) for identifier, _, sequence in read_fasta(fasta)]
        old_records_by_species[species_file] = records
        all_old_ids.extend(identifier for identifier, _ in records)
    if len(all_old_ids) != 1791 or len(set(all_old_ids)) != 1791:
        raise RuntimeError(
            f"Old FASTA IDs expected 1791 unique records, observed "
            f"{len(all_old_ids)} records and {len(set(all_old_ids))} unique IDs"
        )
    absent_new_ids = old_new_ids - set(all_old_ids)
    if absent_new_ids:
        raise RuntimeError(f"Mapped new-predicted IDs absent from FASTA: {sorted(absent_new_ids)}")

    final_protein_records = {
        identifier: (header, sequence.rstrip("*"))
        for identifier, header, sequence in read_fasta(args.retained_proteins)
    }
    deliverable_ids = {row["deliverable_model_id"] for row in retained_rows}
    if set(final_protein_records) != deliverable_ids:
        raise RuntimeError(
            "Retained-protein FASTA and audit deliverable IDs differ: "
            f"FASTA-only={sorted(set(final_protein_records) - deliverable_ids)}; "
            f"audit-only={sorted(deliverable_ids - set(final_protein_records))}"
        )
    invalid = {
        identifier: sorted(set(sequence.upper()) - set("ABCDEFGHIKLMNPQRSTVWXYZOUJB"))
        for identifier, (_, sequence) in final_protein_records.items()
        if set(sequence.upper()) - set("ABCDEFGHIKLMNPQRSTVWXYZOUJB")
    }
    if invalid:
        raise RuntimeError(f"Unexpected amino-acid characters: {invalid}")

    args.output_root.mkdir(parents=True)
    (args.output_root / "source_snapshot").mkdir()
    shutil.copy2(args.audit_tsv, args.output_root / "source_snapshot" / args.audit_tsv.name)
    shutil.copy2(
        args.retained_proteins,
        args.output_root / "source_snapshot" / args.retained_proteins.name,
    )
    shutil.copy2(old_mapping_path, args.output_root / "source_snapshot" / old_mapping_path.name)

    retained_by_species: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in retained_rows:
        retained_by_species[safe_species(row["species"])].append(row)

    manifest_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    expected_totals = {
        "revised_audited_25": 1764,
        "sensitivity_no_independent_new": 1743,
        "sensitivity_annotated_only": 1739,
    }

    for dataset, mode in DATASETS.items():
        outdir = args.output_root / "inputs" / dataset
        outdir.mkdir(parents=True)
        dataset_ids: list[str] = []

        for species_file, old_records in sorted(old_records_by_species.items()):
            annotated = [record for record in old_records if record[0] not in old_new_ids]
            selected_rows: list[dict[str, str]] = []
            if mode == "all_retained":
                selected_rows = retained_by_species.get(species_file, [])
            elif mode == "annotation_corrections_only":
                selected_rows = [
                    row
                    for row in retained_by_species.get(species_file, [])
                    if row["independent_new_gene"] == "No"
                ]

            additions: list[tuple[str, str]] = []
            for row in sorted(selected_rows, key=lambda item: item["deliverable_model_id"]):
                identifier = row["deliverable_model_id"]
                sequence = final_protein_records[identifier][1]
                additions.append((identifier, sequence))
                manifest_rows.append(
                    {
                        "dataset": dataset,
                        "species": row["species"],
                        "record_id": identifier,
                        "source_candidate_id": row["candidate_gene_id"],
                        "record_source": "audited_final_model",
                        "counting_category": row["final_counting_category"],
                        "independent_new_gene": row["independent_new_gene"],
                        "action_relative_to_old_input": "add_final_model",
                        "protein_length_aa": len(sequence),
                    }
                )

            combined = annotated + additions
            write_fasta(outdir / f"{species_file}.fa", combined)
            dataset_ids.extend(identifier for identifier, _ in combined)
            summary_rows.append(
                {
                    "dataset": dataset,
                    "species": species_file,
                    "old_total": len(old_records),
                    "old_candidates_removed": len(old_records) - len(annotated),
                    "annotated_retained": len(annotated),
                    "audited_models_added": len(additions),
                    "final_total": len(combined),
                }
            )

        if len(dataset_ids) != expected_totals[dataset]:
            raise RuntimeError(
                f"{dataset}: expected {expected_totals[dataset]} records, found {len(dataset_ids)}"
            )
        if len(set(dataset_ids)) != len(dataset_ids):
            duplicates = [item for item, count in Counter(dataset_ids).items() if count > 1]
            raise RuntimeError(f"{dataset}: duplicate record IDs: {duplicates}")

    with (args.output_root / "input_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    manifest_fields = [
        "dataset",
        "species",
        "record_id",
        "source_candidate_id",
        "record_source",
        "counting_category",
        "independent_new_gene",
        "action_relative_to_old_input",
        "protein_length_aa",
    ]
    with (args.output_root / "audited_model_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=manifest_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(manifest_rows)

    removal_manifest = []
    audit_by_key = {(row["species"], row["candidate_gene_id"]): row for row in audit_rows}
    for row in old_new_rows:
        audit = audit_by_key[(row["Species"], row["Original_ID"])]
        removal_manifest.append(
            {
                "species": row["Species"],
                "old_orthofinder_id": row["OrthoFinder_ID"],
                "candidate_id": row["Original_ID"],
                "audit_decision": audit["final_counting_category"],
                "retained_in_revised_set": audit["retained_in_revised_bZIP_set"],
                "replacement_model_id": audit["deliverable_model_id"] or "NA",
            }
        )
    with (args.output_root / "old_candidate_removal_and_replacement_manifest.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(removal_manifest[0]), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(removal_manifest)

    files_to_hash = sorted(
        path
        for path in args.output_root.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.tsv"
    )
    with (args.output_root / "SHA256SUMS.tsv").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        handle.write("sha256\trelative_path\n")
        for path in files_to_hash:
            handle.write(f"{sha256(path)}\t{path.relative_to(args.output_root)}\n")

    report = {
        "old_total": 1791,
        "old_annotated": 1739,
        "old_candidates": 52,
        "audit_records": len(audit_rows),
        "retained_final_models": len(retained_rows),
        "independent_new_loci": len(independent),
        "existing_locus_annotation_corrections": len(corrections),
        "excluded_candidate_records": len(audit_rows) - len(retained_rows),
        "dataset_totals": expected_totals,
    }
    (args.output_root / "preparation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
