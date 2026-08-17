#!/usr/bin/env python3
"""Repeat duplication-mode tests after combining cloud with shell for comparability."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


ORDER = ["WGD/segmental", "dispersed", "proximal", "tandem", "singleton"]
SIMULATIONS = 1_000_000
SEED = 20260721


def exact_p(table: np.ndarray) -> tuple[float, int]:
    distribution = stats.random_table(
        table.sum(axis=1), table.sum(axis=0), seed=np.random.default_rng(SEED)
    )
    observed = float(distribution.logpmf(table))
    extreme = 0
    completed = 0
    while completed < SIMULATIONS:
        current = min(20_000, SIMULATIONS - completed)
        sampled = distribution.rvs(size=current)
        extreme += int(np.count_nonzero(distribution.logpmf(sampled) <= observed + 1e-12))
        completed += current
    return (extreme + 1) / (SIMULATIONS + 1), extreme


def metrics(table: np.ndarray) -> dict[str, float | int]:
    chi2, p_value, dof, expected = stats.chi2_contingency(table, correction=False)
    v = math.sqrt(chi2 / (table.sum() * min(table.shape[0] - 1, table.shape[1] - 1)))
    return {
        "chi_square": float(chi2),
        "df": int(dof),
        "pearson_p_value": float(p_value),
        "cramers_v": float(v),
        "minimum_expected_count": float(expected.min()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--analysis-dir", required=True, type=Path)
    args = parser.parse_args()

    reports = {}
    for label in ("fixed_annotated_only", "de_novo_annotated_only"):
        path = args.analysis_dir / f"{label}_duplication_gene_assignments.tsv"
        data = pd.read_csv(path, sep="\t")
        data["Comparable_category"] = data["Category"].replace({"cloud": "shell"})
        full = (
            pd.crosstab(data["Comparable_category"], data["Duplication_type"])
            .reindex(index=["core", "softcore", "shell"], columns=ORDER, fill_value=0)
            .astype(int)
        )
        full.to_csv(args.analysis_dir / f"{label}_duplication_3x5_cloud_with_shell.tsv", sep="\t")
        p_value, extreme = exact_p(full.to_numpy())

        collapsed = full.copy()
        collapsed["local/other"] = collapsed[["proximal", "tandem", "singleton"]].sum(axis=1)
        collapsed = collapsed[["WGD/segmental", "dispersed", "local/other"]]
        collapsed.to_csv(
            args.analysis_dir / f"{label}_duplication_3x3_cloud_with_shell.tsv", sep="\t"
        )
        reports[label] = {
            "genes": int(len(data)),
            "full_3x5": full.to_dict(orient="index"),
            "full_3x5_exact": {
                "test": "Fisher-Freeman-Halton Monte Carlo exact test",
                "simulations": SIMULATIONS,
                "seed": SEED,
                "p_value": float(p_value),
                "extreme_simulations": int(extreme),
                **metrics(full.to_numpy()),
            },
            "collapsed_3x3": {
                "test": "Pearson chi-square test",
                **metrics(collapsed.to_numpy()),
            },
        }

    output = args.analysis_dir / "duplication_cloud_collapsed_comparability.json"
    output.write_text(json.dumps(reports, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(reports, indent=2))


if __name__ == "__main__":
    main()
