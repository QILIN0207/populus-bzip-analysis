#!/usr/bin/env python3
"""Validate clean per-species GFF3/CDS/protein deliverables with gffread."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--genome-root", required=True, type=Path)
    parser.add_argument("--temp-root", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    return parser.parse_args()


def only_file(directory: Path, pattern: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {pattern} in {directory}, found {matches}")
    return matches[0]


def read_fasta(path: Path, trim_stop: bool = False) -> dict[str, str]:
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
                    sequence = "".join(chunks).upper()
                    records[identifier] = sequence.rstrip("*") if trim_stop else sequence
                identifier = line[1:].split()[0]
                chunks = []
            else:
                if identifier is None:
                    raise RuntimeError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if identifier is not None:
        sequence = "".join(chunks).upper()
        records[identifier] = sequence.rstrip("*") if trim_stop else sequence
    return records


def main() -> None:
    args = parse_args()
    if args.output_json.exists():
        raise FileExistsError(f"Refusing to overwrite report: {args.output_json}")
    species_root = args.package_root / "by_species"
    if not species_root.is_dir():
        raise FileNotFoundError(species_root)
    args.temp_root.mkdir(parents=True, exist_ok=True)

    per_species: list[dict[str, object]] = []
    for species_dir in sorted(path for path in species_root.iterdir() if path.is_dir()):
        species = species_dir.name
        gff = only_file(species_dir, "*.gff3")
        expected_cds_path = only_file(species_dir, "*.cds.fasta")
        expected_protein_path = only_file(species_dir, "*.protein.fasta")
        genome = args.genome_root / species / f"{species}.genome.fa"
        if not genome.is_file():
            raise FileNotFoundError(genome)
        output_dir = args.temp_root / species
        output_dir.mkdir(exist_ok=True)
        observed_cds_path = output_dir / "gffread.cds.fasta"
        observed_protein_path = output_dir / "gffread.protein.fasta"
        process = subprocess.run(
            [
                "gffread",
                "-E",
                str(gff),
                "-g",
                str(genome),
                "-x",
                str(observed_cds_path),
                "-y",
                str(observed_protein_path),
            ],
            text=True,
            capture_output=True,
        )
        row: dict[str, object] = {
            "species": species,
            "returncode": process.returncode,
            "stderr": process.stderr.strip(),
        }
        if process.returncode == 0:
            expected_cds = read_fasta(expected_cds_path)
            expected_protein = read_fasta(expected_protein_path, trim_stop=True)
            observed_cds = read_fasta(observed_cds_path)
            observed_protein = read_fasta(observed_protein_path, trim_stop=True)
            cds_bad = sorted(
                identifier
                for identifier in set(expected_cds) & set(observed_cds)
                if expected_cds[identifier] != observed_cds[identifier]
            )
            protein_bad = sorted(
                identifier
                for identifier in set(expected_protein) & set(observed_protein)
                if expected_protein[identifier] != observed_protein[identifier]
            )
            id_sets_equal = (
                set(expected_cds)
                == set(expected_protein)
                == set(observed_cds)
                == set(observed_protein)
            )
            row.update(
                {
                    "models": len(expected_protein),
                    "id_sets_equal": id_sets_equal,
                    "cds_exact": not cds_bad,
                    "protein_exact": not protein_bad,
                    "cds_mismatch_models": cds_bad,
                    "protein_mismatch_models": protein_bad,
                }
            )
        per_species.append(row)

    summary = {
        "species_checked": len(per_species),
        "models_checked": sum(int(row.get("models", 0)) for row in per_species),
        "all_passed": all(
            row["returncode"] == 0
            and bool(row.get("id_sets_equal"))
            and bool(row.get("cds_exact"))
            and bool(row.get("protein_exact"))
            for row in per_species
        ),
    }
    report = {"summary": summary, "per_species": per_species}
    args.output_json.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    if not summary["all_passed"]:
        raise RuntimeError("gffread validation failed")


if __name__ == "__main__":
    main()
