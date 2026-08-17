#!/usr/bin/env python3
"""Create collision-free ParaAT inputs while preserving original identifiers."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(parts)
                identifier = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if identifier is not None:
        records[identifier] = "".join(parts)
    return records


def write_fasta(records: dict[str, str], aliases: dict[str, str], path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for original_id, sequence in records.items():
            handle.write(f">{aliases[original_id]}\n")
            for start in range(0, len(sequence), 60):
                handle.write(sequence[start:start + 60] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--protein", type=Path, required=True)
    parser.add_argument("--cds", type=Path, required=True)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    protein = read_fasta(args.protein)
    cds = read_fasta(args.cds)
    if set(protein) != set(cds):
        raise SystemExit("Protein and CDS identifiers differ")
    original_ids = list(protein)
    aliases = {identifier: f"BZIP{index:06d}X" for index, identifier in enumerate(original_ids, start=1)}
    pair_rows: list[tuple[str, str]] = []
    with args.pairs.open(encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row:
                continue
            if len(row) != 2 or any(identifier not in aliases for identifier in row):
                raise SystemExit(f"Invalid pair row: {row}")
            pair_rows.append((row[0], row[1]))
    write_fasta(protein, aliases, args.output_dir / "paraat_safe.protein.fasta")
    write_fasta(cds, aliases, args.output_dir / "paraat_safe.cds.fasta")
    with (args.output_dir / "paraat_safe.pairs.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        for gene1, gene2 in pair_rows:
            handle.write(f"{aliases[gene1]}\t{aliases[gene2]}\n")
    with (args.output_dir / "paraat_alias_map.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Original_ID", "ParaAT_safe_ID"])
        writer.writerows((identifier, aliases[identifier]) for identifier in original_ids)
    summary = {"sequence_ids": len(original_ids), "unique_safe_ids": len(set(aliases.values())), "pairs": len(pair_rows)}
    (args.output_dir / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
