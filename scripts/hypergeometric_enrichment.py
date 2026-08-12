#!/usr/bin/env python3
"""Run one-sided hypergeometric enrichment tests with BH correction."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from scipy.stats import hypergeom


REQUIRED_COLUMNS = {
    "database",
    "gene_count",
    "foreground_gene_count",
    "term_background_count",
    "background_gene_count",
}


def bh_adjust(values: pd.Series) -> pd.Series:
    order = values.sort_values().index
    ranks = pd.Series(range(1, len(values) + 1), index=order, dtype=float)
    adjusted = (values.loc[order] * len(values) / ranks).iloc[::-1].cummin().iloc[::-1]
    return adjusted.clip(upper=1).reindex(values.index)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    table = pd.read_csv(args.input, sep="\t")
    missing = REQUIRED_COLUMNS - set(table.columns)
    if missing:
        raise ValueError(f"Input is missing columns: {sorted(missing)}")

    table["p_value"] = hypergeom.sf(
        table["gene_count"] - 1,
        table["background_gene_count"],
        table["term_background_count"],
        table["foreground_gene_count"],
    )
    table["adjusted_p_value"] = table.groupby("database", group_keys=False)[
        "p_value"
    ].apply(bh_adjust)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
