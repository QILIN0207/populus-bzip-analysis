#!/usr/bin/env python3
"""Align MCScanX gene-level IDs to the OrthoFinder transcript-level mapping."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", required=True, type=Path)
    parser.add_argument("--duplication", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    mapping = [row for row in read_tsv(args.mapping) if row["Type"] == "annotated"]
    duplication = read_tsv(args.duplication)
    duplication_by_gene = {(row["species"], row["geneID"]): row for row in duplication}

    output_rows = []
    missing = []
    for row in mapping:
        source = duplication_by_gene.get((row["Species"], row["Gene_key"]))
        if source is None:
            missing.append((row["Species"], row["Original_ID"], row["Gene_key"]))
            continue
        revised = dict(source)
        revised["geneID"] = row["Original_ID"]
        output_rows.append(revised)

    # One annotated input was unassigned in the original OrthoFinder result and
    # therefore has no row in the 1,790-gene duplication table.
    if len(output_rows) != 1738 or len(missing) != 1:
        raise RuntimeError(
            f"Expected 1738 mapped rows and one unassigned record; "
            f"observed {len(output_rows)} and {len(missing)}; missing={missing[:5]}"
        )
    keys = [(row["species"], row["geneID"]) for row in output_rows]
    if len(keys) != len(set(keys)):
        raise RuntimeError("Expanded duplication aliases are not unique")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"rows={len(output_rows)} missing_original_unassigned={missing[0]}")


if __name__ == "__main__":
    main()
