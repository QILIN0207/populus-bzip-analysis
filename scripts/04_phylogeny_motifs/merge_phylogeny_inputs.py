#!/usr/bin/env python3
"""Merge revised Populus OGG representatives with fixed Arabidopsis/rice references."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    identifier: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records.append((identifier, "".join(chunks).upper()))
                identifier = line[1:].split()[0]
                chunks = []
            else:
                if identifier is None:
                    raise RuntimeError(f"Sequence before header in {path}")
                chunks.append(line)
    if identifier is not None:
        records.append((identifier, "".join(chunks).upper()))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representatives-fasta", required=True, type=Path)
    parser.add_argument("--representatives-tsv", required=True, type=Path)
    parser.add_argument("--references-fasta", required=True, type=Path)
    parser.add_argument("--merged-fasta", required=True, type=Path)
    parser.add_argument("--source-tsv", required=True, type=Path)
    args = parser.parse_args()
    for output in (args.merged_fasta, args.source_tsv):
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite: {output}")

    with args.representatives_tsv.open(encoding="utf-8", newline="") as handle:
        representative_rows = list(csv.DictReader(handle, delimiter="\t"))
    metadata = {row["Representative_header"]: row for row in representative_rows}
    pop_records = read_fasta(args.representatives_fasta)
    ref_records = read_fasta(args.references_fasta)
    if len(pop_records) != 79:
        raise RuntimeError(f"Expected 79 Populus representatives, found {len(pop_records)}")
    if len(ref_records) != 164:
        raise RuntimeError(f"Expected 164 At/Os reference proteins, found {len(ref_records)}")

    all_records = ref_records + pop_records
    identifiers = [identifier for identifier, _ in all_records]
    if len(identifiers) != len(set(identifiers)):
        raise RuntimeError("Duplicate sequence IDs in merged phylogeny input")
    if any(not sequence for _, sequence in all_records):
        raise RuntimeError("Empty protein sequence in merged phylogeny input")

    args.merged_fasta.parent.mkdir(parents=True, exist_ok=True)
    with args.merged_fasta.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in all_records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start : start + 60] + "\n")

    fields = [
        "Sequence_ID",
        "Source",
        "Subfamily_or_category",
        "OGG",
        "Species",
        "Original_gene",
        "Length",
    ]
    with args.source_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for identifier, sequence in ref_records:
            if identifier.startswith("At"):
                source = "Arabidopsis"
            elif identifier.startswith("Os"):
                source = "Rice"
            else:
                raise RuntimeError(f"Unexpected reference identifier: {identifier}")
            group = identifier.rsplit("_", 1)[-1] if "_" in identifier else "NA"
            writer.writerow(
                {
                    "Sequence_ID": identifier,
                    "Source": source,
                    "Subfamily_or_category": group,
                    "OGG": "NA",
                    "Species": source,
                    "Original_gene": identifier.rsplit("_", 1)[0],
                    "Length": len(sequence),
                }
            )
        for identifier, sequence in pop_records:
            if identifier not in metadata:
                raise RuntimeError(f"Missing representative metadata: {identifier}")
            row = metadata[identifier]
            writer.writerow(
                {
                    "Sequence_ID": identifier,
                    "Source": "Populus",
                    "Subfamily_or_category": row["Category"],
                    "OGG": row["New_Orthogroup"],
                    "Species": row["Species"],
                    "Original_gene": row["Original_gene"],
                    "Length": len(sequence),
                }
            )
    print(f"reference_sequences\t{len(ref_records)}")
    print(f"populus_representatives\t{len(pop_records)}")
    print(f"merged_sequences\t{len(all_records)}")


if __name__ == "__main__":
    main()
