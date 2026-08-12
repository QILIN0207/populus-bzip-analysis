#!/usr/bin/env python3
"""Summarize OGG occupancy and copy-number variation across 19 genomes."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


def members(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def pangenome_class(occupancy: int) -> str:
    if occupancy == 19:
        return "core"
    if occupancy in {17, 18}:
        return "softcore"
    if 2 <= occupancy <= 16:
        return "shell"
    if occupancy == 1:
        return "cloud"
    raise ValueError(f"Invalid occupancy: {occupancy}")


def count_unassigned(path: Path, species: list[str]) -> int:
    table = pd.read_csv(path, sep="\t", dtype=str)
    return sum(len(members(value)) for name in species for value in table[name])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orthogroups", required=True, type=Path)
    parser.add_argument("--unassigned", required=True, type=Path)
    parser.add_argument("--output-classification", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.orthogroups, sep="\t", dtype=str)
    species = [column for column in table.columns if column != "Orthogroup"]
    if len(species) != 19:
        raise ValueError(f"Expected 19 species columns, found {len(species)}")

    rows: list[dict[str, object]] = []
    for record in table.to_dict("records"):
        counts = [len(members(record[name])) for name in species]
        occupancy = sum(count > 0 for count in counts)
        rows.append(
            {
                "orthogroup": record["Orthogroup"],
                "pangenome_class": pangenome_class(occupancy),
                "occupied_genomes": occupancy,
                "total_genes": sum(counts),
                "copy_number_variable": len(set(counts)) > 1,
                **dict(zip(species, counts)),
            }
        )

    classification = pd.DataFrame(rows).sort_values("orthogroup")
    unassigned = count_unassigned(args.unassigned, species)
    class_counts = classification["pangenome_class"].value_counts()
    summary = {
        "input_proteins": int(classification["total_genes"].sum()) + unassigned,
        "assigned_proteins": int(classification["total_genes"].sum()),
        "unassigned_proteins": unassigned,
        "orthogroups": int(len(classification)),
        "core": int(class_counts.get("core", 0)),
        "softcore": int(class_counts.get("softcore", 0)),
        "shell": int(class_counts.get("shell", 0)),
        "cloud": int(class_counts.get("cloud", 0)),
        "copy_number_variable_orthogroups": int(
            classification["copy_number_variable"].sum()
        ),
    }

    args.output_classification.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    classification.to_csv(args.output_classification, sep="\t", index=False)
    args.output_summary.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
