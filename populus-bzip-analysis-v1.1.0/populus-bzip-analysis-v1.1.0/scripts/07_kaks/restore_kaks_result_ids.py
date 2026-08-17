#!/usr/bin/env python3
"""Map collision-free ParaAT/KaKs outputs back to complete original IDs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_aliases(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        return {row["Original_ID"]: row["ParaAT_safe_ID"] for row in rows}


def read_pairs(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--alias-map", type=Path, required=True)
    parser.add_argument("--alias-result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    aliases = read_aliases(args.alias_map)
    pairs = read_pairs(args.pairs)
    missing: list[str] = []
    mapped = 0
    for pair in pairs:
        gene1, gene2 = pair["Gene1"], pair["Gene2"]
        alias_name = f"{aliases[gene1]}-{aliases[gene2]}.cds_aln.kaks.tsv"
        source = args.alias_result_dir / alias_name
        if not source.exists():
            missing.append(alias_name)
            continue
        target = args.output_dir / f"{gene1}-{gene2}.cds_aln.kaks.tsv"
        lines = source.read_text(encoding="utf-8", errors="replace").splitlines()
        if len(lines) != 2:
            raise SystemExit(f"Expected header plus one data row in {source}")
        fields = lines[1].split("\t")
        fields[0] = f"{gene1}-{gene2}"
        target.write_text(lines[0] + "\n" + "\t".join(fields) + "\n", encoding="utf-8")
        mapped += 1
    summary = {"expected": len(pairs), "mapped": mapped, "missing": len(missing)}
    (args.output_dir / "MAPPING_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    if missing:
        raise SystemExit(f"Missing {len(missing)} alias results; first: {missing[:3]}")


if __name__ == "__main__":
    main()
