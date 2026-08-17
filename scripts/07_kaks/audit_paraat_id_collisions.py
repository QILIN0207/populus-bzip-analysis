#!/usr/bin/env python3
"""Audit ParaAT 2.0 ID-version stripping collisions in paired FASTA inputs."""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


VERSION_RE = re.compile(r"\.\d+$")


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts)))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def normalized(identifier: str) -> str:
    return VERSION_RE.sub("", identifier)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    protein_records = parse_fasta(args.protein)
    cds_records = parse_fasta(args.cds)
    protein = dict(protein_records)
    cds = dict(cds_records)
    if set(protein) != set(cds):
        raise SystemExit("Protein and CDS identifiers differ")

    groups: dict[str, list[str]] = defaultdict(list)
    for identifier, _sequence in protein_records:
        groups[normalized(identifier)].append(identifier)
    collisions = {key: ids for key, ids in groups.items() if len(ids) > 1}

    pair_rows: list[tuple[str, str]] = []
    with args.pairs.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            if len(row) != 2:
                raise SystemExit(f"Expected two pair columns, found {len(row)}")
            pair_rows.append((row[0], row[1]))

    used_ids = {identifier for pair in pair_rows for identifier in pair}
    affected_pair_indices: set[int] = set()
    detail_rows: list[dict[str, object]] = []
    for key, ids in sorted(collisions.items()):
        protein_unique = len({protein[identifier] for identifier in ids})
        cds_unique = len({cds[identifier] for identifier in ids})
        used = [identifier for identifier in ids if identifier in used_ids]
        for index, pair in enumerate(pair_rows, start=1):
            if any(identifier in pair for identifier in ids):
                affected_pair_indices.add(index)
        detail_rows.append(
            {
                "normalized_id": key,
                "original_ids": ";".join(ids),
                "n_original_ids": len(ids),
                "n_protein_sequences": protein_unique,
                "n_cds_sequences": cds_unique,
                "used_ids": ";".join(used),
                "n_used_ids": len(used),
                "sequences_identical": protein_unique == 1 and cds_unique == 1,
            }
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    detail_path = args.output_dir / "paraat_version_collision_details.tsv"
    with detail_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(detail_rows[0]) if detail_rows else ["normalized_id"], delimiter="\t")
        writer.writeheader()
        writer.writerows(detail_rows)

    summary = {
        "protein_records": len(protein_records),
        "cds_records": len(cds_records),
        "paraat_normalized_unique_ids": len(groups),
        "collision_groups": len(collisions),
        "original_ids_in_collision_groups": sum(len(ids) for ids in collisions.values()),
        "collision_groups_with_nonidentical_sequences": sum(
            not bool(row["sequences_identical"]) for row in detail_rows
        ),
        "pairs_total": len(pair_rows),
        "pairs_touching_collision_groups": len(affected_pair_indices),
    }
    (args.output_dir / "paraat_version_collision_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
