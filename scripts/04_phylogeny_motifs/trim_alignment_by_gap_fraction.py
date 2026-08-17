#!/usr/bin/env python3
"""Remove alignment columns whose gap fraction is greater than a threshold."""

from __future__ import annotations

import argparse
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
                    records.append((identifier, "".join(chunks)))
                identifier = line[1:].split()[0]
                chunks = []
            else:
                if identifier is None:
                    raise RuntimeError(f"Sequence before header in {path}")
                chunks.append(line)
    if identifier is not None:
        records.append((identifier, "".join(chunks)))
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--alignment", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--max-gap-fraction", type=float, default=0.80)
    args = parser.parse_args()
    for output in (args.output, args.report):
        if output.exists():
            raise FileExistsError(f"Refusing to overwrite: {output}")
    if not 0 <= args.max_gap_fraction <= 1:
        raise ValueError("--max-gap-fraction must be between 0 and 1")

    records = read_fasta(args.alignment)
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise RuntimeError(f"Alignment lengths differ: {sorted(lengths)}")
    alignment_length = lengths.pop()
    sequence_count = len(records)
    keep_columns = []
    for column in range(alignment_length):
        gap_count = sum(sequence[column] in "-." for _, sequence in records)
        if gap_count / sequence_count <= args.max_gap_fraction:
            keep_columns.append(column)

    with args.output.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            trimmed = "".join(sequence[index] for index in keep_columns)
            handle.write(f">{identifier}\n")
            for start in range(0, len(trimmed), 60):
                handle.write(trimmed[start : start + 60] + "\n")

    retained = len(keep_columns)
    removed = alignment_length - retained
    with args.report.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("Metric\tValue\n")
        handle.write(f"Sequence_number\t{sequence_count}\n")
        handle.write(f"Original_alignment_columns\t{alignment_length}\n")
        handle.write(f"Gap_fraction_removal_rule\t>{args.max_gap_fraction:.2f}\n")
        handle.write(f"Retained_columns_gap_le_{args.max_gap_fraction:.2f}\t{retained}\n")
        handle.write(f"Removed_columns_gap_gt_{args.max_gap_fraction:.2f}\t{removed}\n")
        handle.write(f"Removed_column_percent\t{removed / alignment_length * 100:.2f}\n")
    print(f"sequence_number\t{sequence_count}")
    print(f"original_columns\t{alignment_length}")
    print(f"retained_columns\t{retained}")
    print(f"removed_columns\t{removed}")


if __name__ == "__main__":
    main()
