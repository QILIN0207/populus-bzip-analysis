#!/usr/bin/env python3
"""Assign revised PtbZIP names to the 79 primary de novo OGGs."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


CATEGORY_ORDER = ["core", "softcore", "shell", "cloud"]
CATEGORY_PREFIX = {
    "core": "CR",
    "softcore": "SC",
    "shell": "SH",
    "cloud": "CL",
}


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
                    records[identifier] = "".join(chunks).upper()
                identifier = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if identifier is not None:
        records[identifier] = "".join(chunks).upper()
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")


def ogg_number(ogg: str) -> int:
    if not ogg.startswith("OG"):
        raise ValueError(f"Unexpected OGG ID: {ogg}")
    return int(ogg[2:])


def split_values(value: str) -> list[str]:
    return [item for item in value.split(";") if item and item != "NA"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--representatives-tsv", required=True, type=Path)
    parser.add_argument("--representatives-fasta", required=True, type=Path)
    parser.add_argument("--old-to-new-mapping", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    classification = read_tsv(args.classification)
    if len(classification) != 79:
        raise RuntimeError(f"Expected 79 revised OGGs, found {len(classification)}")
    counts = Counter(row["Category"] for row in classification)
    expected_counts = {"core": 43, "softcore": 20, "shell": 15, "cloud": 1}
    if dict(counts) != expected_counts:
        raise RuntimeError(f"Unexpected revised category counts: {dict(counts)}")

    ordered_rows: list[dict[str, str]] = []
    for category in CATEGORY_ORDER:
        ordered_rows.extend(
            sorted(
                (row for row in classification if row["Category"] == category),
                key=lambda row: ogg_number(row["Orthogroup"]),
            )
        )

    nomenclature_rows: list[dict[str, object]] = []
    pangene_by_ogg: dict[str, str] = {}
    for ordinal, row in enumerate(ordered_rows, start=1):
        category = row["Category"]
        pangene = f"PtbZIP.{CATEGORY_PREFIX[category]}{ordinal:03d}"
        ogg = row["Orthogroup"]
        pangene_by_ogg[ogg] = pangene
        nomenclature_rows.append(
            {
                "Pangene_ID": pangene,
                "New_Orthogroup": ogg,
                "Category": category,
                "Species_present": row["Species_present"],
                "Total_genes": row["Total_genes"],
                "Copy_number_variable": row["Copy_number_variable"],
                "Old_orthogroup_contributors": row.get(
                    "Old_orthogroup_contributors", "NA"
                ),
                "Old_pangene_contributors": row.get("Old_pangene_contributors", "NA"),
                "Merge_of_multiple_old_OGGs": row.get(
                    "Merge_of_multiple_old_OGGs", "NA"
                ),
                "Naming_rule": (
                    "category_order_core_softcore_shell_cloud_then_new_OGG_numeric_order"
                ),
            }
        )

    expected_ranges = {
        "core": ("PtbZIP.CR001", "PtbZIP.CR043"),
        "softcore": ("PtbZIP.SC044", "PtbZIP.SC063"),
        "shell": ("PtbZIP.SH064", "PtbZIP.SH078"),
        "cloud": ("PtbZIP.CL079", "PtbZIP.CL079"),
    }
    for category, (first_expected, last_expected) in expected_ranges.items():
        names = [
            str(row["Pangene_ID"])
            for row in nomenclature_rows
            if row["Category"] == category
        ]
        if (names[0], names[-1]) != (first_expected, last_expected):
            raise RuntimeError(f"Unexpected {category} name range: {names[0]}-{names[-1]}")

    representative_rows = read_tsv(args.representatives_tsv)
    representative_sequences = read_fasta(args.representatives_fasta)
    if len(representative_rows) != 79 or len(representative_sequences) != 79:
        raise RuntimeError("Expected 79 representative rows and sequences")
    representative_by_ogg = {
        row["New_Orthogroup"]: row for row in representative_rows
    }
    pangene_representative_rows: list[dict[str, object]] = []
    pangene_fasta: list[tuple[str, str]] = []
    for nomenclature in nomenclature_rows:
        ogg = str(nomenclature["New_Orthogroup"])
        pangene = str(nomenclature["Pangene_ID"])
        representative = dict(representative_by_ogg[ogg])
        old_header = representative["Representative_header"]
        if old_header not in representative_sequences:
            raise RuntimeError(f"Missing representative sequence: {old_header}")
        suffix = old_header.split("__", 1)[1]
        new_header = f"{pangene}__{suffix}"
        pangene_representative_rows.append(
            {
                "Pangene_ID": pangene,
                **representative,
                "Representative_header": new_header,
            }
        )
        pangene_fasta.append((new_header, representative_sequences[old_header]))

    mapping_rows = [
        row
        for row in read_tsv(args.old_to_new_mapping)
        if row.get("dataset") == "revised_audited_25"
    ]
    if len(mapping_rows) != 86:
        raise RuntimeError(f"Expected 86 historical mapping rows, found {len(mapping_rows)}")
    old_to_new_rows: list[dict[str, object]] = []
    for row in mapping_rows:
        new_oggs = split_values(row.get("Mapped_new_OGGs", ""))
        new_names = [pangene_by_ogg[ogg] for ogg in new_oggs]
        old_to_new_rows.append(
            {
                "Old_Orthogroup": row["Old_Orthogroup"],
                "Old_PtbZIP_ID": row["Old_pangeneID"],
                "Old_category": row["Old_category"],
                "Surviving_annotated_members": row["Surviving_annotated_members"],
                "Mapped_new_OGGs": ";".join(new_oggs) or "NA",
                "New_PtbZIP_IDs": ";".join(new_names) or "NA",
                "Mapping_status": row["Mapping_status"],
            }
        )

    write_tsv(
        args.output_dir / "New79_PtbZIP_nomenclature.tsv",
        nomenclature_rows,
        list(nomenclature_rows[0]),
    )
    write_tsv(
        args.output_dir / "New79_PtbZIP_representatives.tsv",
        pangene_representative_rows,
        list(pangene_representative_rows[0]),
    )
    write_fasta(
        args.output_dir / "New79_PtbZIP_representatives.full_length.fasta",
        pangene_fasta,
    )
    write_tsv(
        args.output_dir / "Old86_to_New79_PtbZIP_nomenclature_mapping.tsv",
        old_to_new_rows,
        list(old_to_new_rows[0]),
    )

    summary = {
        "total_revised_pangenes": len(nomenclature_rows),
        "category_counts": expected_counts,
        "name_ranges": expected_ranges,
        "naming_order": CATEGORY_ORDER,
        "within_category_order": "ascending revised Orthogroup numeric ID",
        "unique_pangene_names": len({row["Pangene_ID"] for row in nomenclature_rows}),
        "representatives_renamed": len(pangene_fasta),
        "old_to_new_mapping_rows": len(old_to_new_rows),
        "all_checks_passed": True,
    }
    (args.output_dir / "nomenclature_summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
