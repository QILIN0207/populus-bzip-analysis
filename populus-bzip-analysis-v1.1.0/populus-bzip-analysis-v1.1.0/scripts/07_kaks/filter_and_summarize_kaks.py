#!/usr/bin/env python3
"""Merge, filter, annotate, and statistically summarize current79 Ka/Ks results."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy import stats


RAW_COLUMNS = [
    "Sequence", "Method", "Ka", "Ks", "Ka/Ks", "P-Value(Fisher)", "Length", "S-Sites", "N-Sites",
    "Fold-Sites(0:2:4)", "Substitutions", "S-Substitutions", "N-Substitutions",
    "Fold-S-Substitutions(0:2:4)", "Fold-N-Substitutions(0:2:4)", "Divergence-Time",
    "Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)", "GC(1:2:3)", "ML-Score", "AICc",
    "Akaike-Weight", "Model",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def numeric(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text or text.upper() in {"NA", "NAN", "-NAN", "INF", "-INF"}:
        return None
    try:
        number = float(text)
    except ValueError:
        return None
    return number if math.isfinite(number) else None


def write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def cliffs_delta(group_a: list[float], group_b: list[float]) -> float:
    greater = sum(a > b for a in group_a for b in group_b)
    less = sum(a < b for a in group_a for b in group_b)
    return (greater - less) / (len(group_a) * len(group_b))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit(f"Refusing to overwrite existing output directory: {args.output_dir}")
    args.output_dir.mkdir(parents=True)

    pairs = read_tsv(args.pairs)
    rows: list[dict[str, object]] = []
    missing_outputs: list[str] = []
    for pair in pairs:
        expected_name = f"{pair['Gene1']}-{pair['Gene2']}.cds_aln.kaks.tsv"
        result_path = args.result_dir / expected_name
        raw: dict[str, str] = {column: "NA" for column in RAW_COLUMNS}
        if result_path.exists():
            result_rows = read_tsv(result_path)
            if len(result_rows) != 1:
                raise ValueError(f"Expected one Ka/Ks row in {result_path}; observed {len(result_rows)}")
            raw.update(result_rows[0])
        else:
            missing_outputs.append(expected_name)

        ka = numeric(raw.get("Ka"))
        ks = numeric(raw.get("Ks"))
        ratio = numeric(raw.get("Ka/Ks"))
        reasons: list[str] = []
        if ka is None or ks is None or ratio is None or ka < 0 or ks <= 0 or ratio < 0:
            reasons.append("undefined_or_nonfinite_estimate")
        else:
            if ks < 0.01:
                reasons.append("Ks_lt_0.01")
            if ks > 2:
                reasons.append("Ks_gt_2")
            if ka > 2:
                reasons.append("Ka_gt_2")
            if ratio > 10:
                reasons.append("KaKs_gt_10")

        retained = not reasons and result_path.exists()
        annotated: dict[str, object] = dict(raw)
        annotated["KaKs"] = raw.get("Ka/Ks", "NA")
        annotated.update({
            "Pair_ID": pair["Pair_ID"],
            "Orthogroup": pair["Orthogroup"],
            "Pangene_ID": pair["Pangene_ID"],
            "Gene1": pair["Gene1"],
            "Gene2": pair["Gene2"],
            "Gene1_species": pair["Gene1_species"],
            "Gene2_species": pair["Gene2_species"],
            "Gene1_source_status": pair["Gene1_source_status"],
            "Gene2_source_status": pair["Gene2_source_status"],
            "Subfamily": pair["Subfamily"],
            "PangenomeClass": pair["Pangenome_class"],
            "CoreStatus": "Core" if pair["Pangenome_class"] == "core" else "Non-core",
            "mapping_version": "New79_20260808",
            "Analysis_status": "retained" if retained else "excluded",
            "Exclusion_reasons": ";".join(reasons) if reasons else "NA",
            "Primary_exclusion_reason": reasons[0] if reasons else "NA",
            "Selection": (
                "Ka/Ks <= 1" if retained and ratio is not None and ratio <= 1
                else "Ka/Ks > 1" if retained and ratio is not None
                else "NA"
            ),
        })
        rows.append(annotated)

    if missing_outputs:
        raise ValueError(f"Missing {len(missing_outputs)} Ka/Ks outputs; first: {missing_outputs[:3]}")

    raw_fieldnames = RAW_COLUMNS + [
        "Pair_ID", "Orthogroup", "Pangene_ID", "Gene1", "Gene2", "Gene1_species", "Gene2_species",
        "Gene1_source_status", "Gene2_source_status", "Subfamily", "PangenomeClass", "CoreStatus",
        "mapping_version", "Analysis_status", "Exclusion_reasons", "Primary_exclusion_reason", "Selection",
    ]
    filtered = [row for row in rows if row["Analysis_status"] == "retained"]
    excluded = [row for row in rows if row["Analysis_status"] == "excluded"]
    write_rows(args.output_dir / "current79_kaks_all_attempted.tsv", rows, raw_fieldnames)
    write_rows(args.output_dir / "current79_kaks_filtered.tsv", filtered, raw_fieldnames)
    write_rows(args.output_dir / "current79_kaks_excluded.tsv", excluded, raw_fieldnames)

    exclusion_primary = Counter(str(row["Primary_exclusion_reason"]) for row in excluded)
    exclusion_any = Counter(
        reason
        for row in excluded
        for reason in str(row["Exclusion_reasons"]).split(";")
        if reason and reason != "NA"
    )
    ratio_values = [float(row["Ka/Ks"]) for row in filtered]
    selection_counts = Counter(str(row["Selection"]) for row in filtered)

    summary_rows: list[dict[str, object]] = [
        {"Metric": "Attempted comparisons", "Value": len(rows), "Notes": "Unique unordered anchor-member pairs"},
        {"Metric": "Retained after filtering", "Value": len(filtered), "Notes": "Used in quantitative summaries and plots"},
        {"Metric": "Excluded comparisons", "Value": len(excluded), "Notes": "Listed separately with explicit reasons"},
        {"Metric": "Retained Ka/Ks <= 1", "Value": selection_counts.get("Ka/Ks <= 1", 0), "Notes": "Purifying/neutral category"},
        {"Metric": "Retained Ka/Ks > 1", "Value": selection_counts.get("Ka/Ks > 1", 0), "Notes": "Elevated-ratio category"},
        {"Metric": "Median retained Ka/Ks", "Value": float(np.median(ratio_values)), "Notes": "Across retained comparisons"},
    ]
    for reason in ("undefined_or_nonfinite_estimate", "Ks_lt_0.01", "Ks_gt_2", "Ka_gt_2", "KaKs_gt_10"):
        summary_rows.append({
            "Metric": f"Primary exclusion: {reason}",
            "Value": exclusion_primary.get(reason, 0),
            "Notes": "Mutually exclusive first applicable reason",
        })
        summary_rows.append({
            "Metric": f"Any exclusion flag: {reason}",
            "Value": exclusion_any.get(reason, 0),
            "Notes": "A comparison may have multiple flags",
        })
    write_rows(args.output_dir / "current79_kaks_filter_summary.tsv", summary_rows, ["Metric", "Value", "Notes"])

    group_summary_rows: list[dict[str, object]] = []
    for field in ("PangenomeClass", "CoreStatus", "Subfamily", "Pangene_ID"):
        grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
        for row in filtered:
            grouped[str(row[field])].append(row)
        for group, group_rows in sorted(grouped.items()):
            values = [float(row["Ka/Ks"]) for row in group_rows]
            group_summary_rows.append({
                "Grouping_variable": field,
                "Group": group,
                "n": len(values),
                "Median_KaKs": float(np.median(values)),
                "Mean_KaKs": float(np.mean(values)),
                "Q1_KaKs": float(np.quantile(values, 0.25)),
                "Q3_KaKs": float(np.quantile(values, 0.75)),
                "KaKs_le_1_n": sum(value <= 1 for value in values),
                "KaKs_gt_1_n": sum(value > 1 for value in values),
            })
    write_rows(
        args.output_dir / "current79_kaks_group_summary.tsv",
        group_summary_rows,
        ["Grouping_variable", "Group", "n", "Median_KaKs", "Mean_KaKs", "Q1_KaKs", "Q3_KaKs", "KaKs_le_1_n", "KaKs_gt_1_n"],
    )

    test_rows: list[dict[str, object]] = []
    core = [float(row["Ka/Ks"]) for row in filtered if row["CoreStatus"] == "Core"]
    noncore = [float(row["Ka/Ks"]) for row in filtered if row["CoreStatus"] == "Non-core"]
    mann = stats.mannwhitneyu(core, noncore, alternative="two-sided", method="asymptotic")
    test_rows.append({
        "Test": "Core versus non-core Ka/Ks",
        "Method": "two-sided Mann-Whitney U",
        "Statistic": float(mann.statistic),
        "P_value": float(mann.pvalue),
        "Effect_size": cliffs_delta(core, noncore),
        "Effect_size_name": "Cliff's delta (Core minus Non-core)",
        "Group_sizes": f"{len(core)};{len(noncore)}",
    })

    class_groups = [
        [float(row["Ka/Ks"]) for row in filtered if row["PangenomeClass"] == category]
        for category in ("core", "softcore", "shell", "cloud")
    ]
    class_groups_nonempty = [group for group in class_groups if group]
    kruskal_class = stats.kruskal(*class_groups_nonempty)
    test_rows.append({
        "Test": "Ka/Ks among pangenome classes",
        "Method": "Kruskal-Wallis",
        "Statistic": float(kruskal_class.statistic),
        "P_value": float(kruskal_class.pvalue),
        "Effect_size": "NA",
        "Effect_size_name": "NA",
        "Group_sizes": ";".join(str(len(group)) for group in class_groups),
    })

    subfamily_groups: dict[str, list[float]] = defaultdict(list)
    for row in filtered:
        subfamily_groups[str(row["Subfamily"])].append(float(row["Ka/Ks"]))
    eligible_subfamilies = {name: values for name, values in subfamily_groups.items() if len(values) >= 2}
    kruskal_subfamily = stats.kruskal(*eligible_subfamilies.values())
    test_rows.append({
        "Test": "Ka/Ks among subfamilies",
        "Method": "Kruskal-Wallis",
        "Statistic": float(kruskal_subfamily.statistic),
        "P_value": float(kruskal_subfamily.pvalue),
        "Effect_size": "NA",
        "Effect_size_name": "NA",
        "Group_sizes": ";".join(f"{name}:{len(values)}" for name, values in sorted(eligible_subfamilies.items())),
    })
    write_rows(
        args.output_dir / "current79_kaks_statistical_tests.tsv",
        test_rows,
        ["Test", "Method", "Statistic", "P_value", "Effect_size", "Effect_size_name", "Group_sizes"],
    )

    summary = {
        "attempted": len(rows),
        "retained": len(filtered),
        "excluded": len(excluded),
        "selection_counts": dict(selection_counts),
        "median_KaKs": float(np.median(ratio_values)),
        "primary_exclusion_counts": dict(exclusion_primary),
        "any_exclusion_flag_counts": dict(exclusion_any),
        "core_n": len(core),
        "noncore_n": len(noncore),
        "tests": test_rows,
    }
    (args.output_dir / "SUMMARY.json").write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=True))


if __name__ == "__main__":
    main()
