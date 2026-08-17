#!/usr/bin/env python3
"""Calculate TE overlap and nonredundant coverage for gene analysis windows."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


GENE_COLUMNS = {
    "gene_id", "species", "chrom", "window_start", "window_end"
}
TE_COLUMNS = {
    "te_id", "species", "chrom", "te_start", "te_end", "te_class", "te_subclass"
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", required=True, type=Path)
    parser.add_argument("--tes", required=True, type=Path)
    parser.add_argument("--output-events", required=True, type=Path)
    parser.add_argument("--output-summary", required=True, type=Path)
    return parser.parse_args()


def require_columns(table: pd.DataFrame, required: set[str], label: str) -> None:
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{label} is missing columns: {sorted(missing)}")


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = 0
    current_start, current_end = sorted(intervals)[0]
    for start, end in sorted(intervals)[1:]:
        if start <= current_end:
            current_end = max(current_end, end)
        else:
            merged += current_end - current_start
            current_start, current_end = start, end
    return merged + current_end - current_start


def main() -> None:
    args = parse_args()
    genes = pd.read_csv(args.genes, sep="\t")
    tes = pd.read_csv(args.tes, sep="\t")
    require_columns(genes, GENE_COLUMNS, "gene table")
    require_columns(tes, TE_COLUMNS, "TE table")

    for start, end, label, table in [
        ("window_start", "window_end", "gene", genes),
        ("te_start", "te_end", "TE", tes),
    ]:
        if (table[start] < 0).any() or (table[end] <= table[start]).any():
            raise ValueError(f"Invalid zero-based half-open {label} interval")

    if genes["gene_id"].duplicated().any():
        raise ValueError("gene_id values must be unique")
    if tes["te_id"].duplicated().any():
        raise ValueError("te_id values must be unique")

    te_groups = {
        key: frame.sort_values("te_start")
        for key, frame in tes.groupby(["species", "chrom"], sort=False)
    }
    event_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for gene in genes.to_dict("records"):
        window_start = int(gene["window_start"])
        window_end = int(gene["window_end"])
        group = te_groups.get((gene["species"], gene["chrom"]))
        covered: list[tuple[int, int]] = []
        event_count = 0

        if group is not None:
            candidates = group[
                (group["te_start"] < window_end) & (group["te_end"] > window_start)
            ]
            for te in candidates.to_dict("records"):
                overlap_start = max(window_start, int(te["te_start"]))
                overlap_end = min(window_end, int(te["te_end"]))
                covered.append((overlap_start, overlap_end))
                event_count += 1
                event_rows.append(
                    {
                        "gene_id": gene["gene_id"],
                        "te_id": te["te_id"],
                        "species": gene["species"],
                        "chrom": gene["chrom"],
                        "overlap_start": overlap_start,
                        "overlap_end": overlap_end,
                        "overlap_bp": overlap_end - overlap_start,
                        "te_class": te["te_class"],
                        "te_subclass": te["te_subclass"],
                    }
                )

        covered_bp = union_length(covered)
        window_bp = window_end - window_start
        summary_rows.append(
            {
                "gene_id": gene["gene_id"],
                "species": gene["species"],
                "chrom": gene["chrom"],
                "window_start": window_start,
                "window_end": window_end,
                "window_bp": window_bp,
                "te_events": event_count,
                "te_positive": event_count > 0,
                "te_covered_bp": covered_bp,
                "te_covered_fraction": covered_bp / window_bp,
            }
        )

    args.output_events.parent.mkdir(parents=True, exist_ok=True)
    args.output_summary.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(event_rows).to_csv(args.output_events, sep="\t", index=False)
    pd.DataFrame(summary_rows).to_csv(args.output_summary, sep="\t", index=False)


if __name__ == "__main__":
    main()
