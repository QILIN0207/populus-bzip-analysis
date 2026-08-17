#!/usr/bin/env python3
"""Merge two-column HTSeq count files into one integer count matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    files = sorted(args.counts_root.glob("*/counts.tsv"))
    if not files:
        raise FileNotFoundError(f"No */counts.tsv files under {args.counts_root}")
    merged: pd.DataFrame | None = None
    for path in files:
        sample = path.parent.name
        frame = pd.read_csv(path, sep="\t", header=None, names=["gene_id", sample])
        frame = frame[~frame["gene_id"].astype(str).str.startswith("__")]
        merged = frame if merged is None else merged.merge(frame, on="gene_id", how="outer")
    assert merged is not None
    merged = merged.fillna(0)
    merged.iloc[:, 1:] = merged.iloc[:, 1:].astype(int)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()

