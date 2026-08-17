#!/usr/bin/env python3
"""Input utilities shared by the transposable-element analyses."""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


EXCLUDED_REPEAT_CLASSES = {
    "Simple_repeat",
    "Low_complexity",
    "Satellite",
    "rRNA",
    "scRNA",
    "snRNA",
    "srpRNA",
    "tRNA",
    "RNA",
}


@dataclass(frozen=True)
class TERecord:
    """A zero-based, half-open RepeatMasker interval."""

    start0: int
    end0: int
    repeat_name: str
    repeat_class: str


def read_seq_mapping(path: Path) -> dict[str, str]:
    """Read original-to-normalized sequence identifiers."""
    result: dict[str, str] = {}
    for line in path.open(encoding="utf-8", errors="replace"):
        parts = line.rstrip("\n").split("\t")
        if not parts:
            continue
        normalized = parts[-1]
        result[parts[0].lstrip(">").split()[0]] = normalized
        for part in parts[1:]:
            if part.startswith("OriSeqID="):
                result[part.split("=", 1)[1]] = normalized
    return result


def read_fai(path: Path) -> dict[str, int]:
    """Read sequence lengths from a FASTA index."""
    result: dict[str, int] = {}
    for line in path.open(encoding="utf-8"):
        parts = line.rstrip("\n").split("\t")
        result[parts[0]] = int(parts[1])
    return result


def read_repeatmasker(path: Path, wanted: set[str]) -> dict[str, list[TERecord]]:
    """Read classified RepeatMasker records for selected sequences."""
    records: dict[str, list[TERecord]] = defaultdict(list)
    for line in path.open(encoding="utf-8", errors="replace"):
        if not re.match(r"^\s*\d+", line):
            continue
        parts = re.split(r"\s+", line.strip())
        if len(parts) < 11 or parts[4] not in wanted:
            continue
        try:
            left, right = int(parts[5]), int(parts[6])
        except ValueError:
            continue
        repeat_class = parts[10]
        if repeat_class in EXCLUDED_REPEAT_CLASSES:
            continue
        start0, end0 = min(left, right) - 1, max(left, right)
        if end0 > start0:
            records[parts[4]].append(
                TERecord(start0, end0, parts[9], repeat_class)
            )
    for chrom in records:
        records[chrom].sort(
            key=lambda record: (record.start0, record.end0, record.repeat_name)
        )
    return records


def merge_coverage(intervals: Iterable[tuple[int, int]]) -> int:
    """Return the covered length after merging overlapping intervals."""
    ordered = sorted(intervals)
    if not ordered:
        return 0
    total = 0
    left, right = ordered[0]
    for start, end in ordered[1:]:
        if start <= right:
            right = max(right, end)
        else:
            total += right - left
            left, right = start, end
    return total + right - left
