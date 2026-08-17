#!/usr/bin/env python3
"""Verify the retained GO/KEGG hypergeometric and BH-adjusted values."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom


def bh_adjust(p_values: list[float]) -> list[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    adjusted = [1.0] * len(p_values)
    running = 1.0
    total = len(p_values)
    for rank_index in range(total - 1, -1, -1):
        original_index, p_value = indexed[rank_index]
        rank = rank_index + 1
        running = min(running, p_value * total / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    table = pd.read_csv(args.input, sep="\t")
    table["recalculated_p_value"] = hypergeom.sf(
        table["gene_count"] - 1,
        table["background_gene_count"],
        table["term_background_count"],
        table["foreground_gene_count"],
    )
    table["recalculated_FDR"] = 1.0
    for _, indices in table.groupby("database", sort=False).groups.items():
        selected = list(indices)
        adjusted = bh_adjust(table.loc[selected, "recalculated_p_value"].astype(float).tolist())
        table.loc[selected, "recalculated_FDR"] = adjusted

    p_diff = (table["recalculated_p_value"] - table["p_value"]).abs()
    fdr_diff = (table["recalculated_FDR"] - table["FDR_adjusted_p_value"]).abs()
    qa = {
        "row_count": int(len(table)),
        "databases": table["database"].value_counts().to_dict(),
        "maximum_p_value_difference": float(p_diff.max()),
        "maximum_FDR_difference": float(fdr_diff.max()),
        "all_p_values_verified": bool((p_diff < 1e-12).all()),
        "all_FDR_values_verified": bool((fdr_diff < 1e-12).all()),
    }
    qa["all_passed"] = qa["all_p_values_verified"] and qa["all_FDR_values_verified"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(qa, indent=2), encoding="utf-8")
    print(json.dumps(qa, indent=2))


if __name__ == "__main__":
    main()
