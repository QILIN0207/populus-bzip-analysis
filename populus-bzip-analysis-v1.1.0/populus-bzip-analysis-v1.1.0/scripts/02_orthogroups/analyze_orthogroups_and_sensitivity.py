#!/usr/bin/env python3
"""Compare audited Populus bZIP OGG runs and reviewer sensitivity analyses."""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


SPECIES_N = 19
CATEGORY_ORDER = ["core", "softcore", "shell", "cloud", "unknown"]
DUPLICATION_ORDER = ["WGD/segmental", "dispersed", "proximal", "tandem", "singleton"]
SIMULATIONS = 1_000_000
SEED = 20260721


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--revision-root", required=True, type=Path)
    parser.add_argument("--old-root", required=True, type=Path)
    parser.add_argument("--audit-tsv", required=True, type=Path)
    parser.add_argument("--duplication-tsv", required=True, type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def split_cell(value: str | float | None) -> list[str]:
    if value is None or pd.isna(value) or not str(value).strip():
        return []
    return [item.strip() for item in str(value).split(",") if item.strip()]


def classify(species_present: int) -> str:
    if species_present == 19:
        return "core"
    if 17 <= species_present <= 18:
        return "softcore"
    if 2 <= species_present <= 16:
        return "shell"
    if species_present == 1:
        return "cloud"
    return "unknown"


def result_dir(run_dir: Path) -> Path:
    candidates = sorted(run_dir.glob("Results_*"))
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one Results_* directory in {run_dir}, found {candidates}")
    return candidates[0]


def load_orthogroups(results: Path) -> tuple[list[str], dict[str, dict[str, list[str]]], list[str]]:
    table = pd.read_csv(results / "Orthogroups" / "Orthogroups.tsv", sep="\t", dtype=str)
    species = list(table.columns[1:])
    if len(species) != SPECIES_N:
        raise RuntimeError(f"Expected 19 species columns, found {len(species)}")
    groups: dict[str, dict[str, list[str]]] = {}
    for _, row in table.iterrows():
        orthogroup = str(row["Orthogroup"])
        groups[orthogroup] = {name: split_cell(row[name]) for name in species}

    unassigned_path = results / "Orthogroups" / "Orthogroups_UnassignedGenes.tsv"
    unassigned: list[str] = []
    if unassigned_path.exists():
        unassigned_table = pd.read_csv(unassigned_path, sep="\t", dtype=str)
        for name in species:
            if name in unassigned_table:
                for value in unassigned_table[name]:
                    unassigned.extend(split_cell(value))
    return species, groups, unassigned


def group_members(group: dict[str, list[str]]) -> set[str]:
    return {identifier for members in group.values() for identifier in members}


def group_stats(
    groups: dict[str, dict[str, list[str]]], species: list[str]
) -> pd.DataFrame:
    rows = []
    for orthogroup, per_species in groups.items():
        counts = [len(per_species.get(name, [])) for name in species]
        present = sum(value > 0 for value in counts)
        rows.append(
            {
                "Orthogroup": orthogroup,
                "Category": classify(present),
                "Species_present": present,
                "Total_genes": sum(counts),
                "Copy_number_variable": len(set(counts)) > 1,
                "Min_copy_number": min(counts),
                "Max_copy_number": max(counts),
                **{name: counts[index] for index, name in enumerate(species)},
            }
        )
    return pd.DataFrame(rows).sort_values("Orthogroup").reset_index(drop=True)


def summary_row(
    dataset: str,
    analysis_mode: str,
    stats_table: pd.DataFrame,
    input_total: int,
    unassigned_count: int,
) -> dict[str, object]:
    row: dict[str, object] = {
        "dataset": dataset,
        "analysis_mode": analysis_mode,
        "input_genes": input_total,
        "assigned_genes": int(stats_table["Total_genes"].sum()),
        "unassigned_genes": unassigned_count,
        "nonempty_OGGs": int(len(stats_table)),
        "copy_number_variable_OGGs": int(stats_table["Copy_number_variable"].sum()),
        "invariant_OGGs": int((~stats_table["Copy_number_variable"]).sum()),
    }
    for category in CATEGORY_ORDER:
        subset = stats_table[stats_table["Category"] == category]
        row[f"{category}_OGGs"] = int(len(subset))
        row[f"{category}_genes"] = int(subset["Total_genes"].sum())
    return row


def copy_groups(groups: dict[str, dict[str, list[str]]]) -> dict[str, dict[str, list[str]]]:
    return {
        orthogroup: {species: list(members) for species, members in per_species.items()}
        for orthogroup, per_species in groups.items()
    }


def fixed_membership_dataset(
    old_groups: dict[str, dict[str, list[str]]],
    old_new_ids: set[str],
    audit_rows: list[dict[str, str]],
    mode: str,
) -> dict[str, dict[str, list[str]]]:
    groups = copy_groups(old_groups)
    for per_species in groups.values():
        for species, members in per_species.items():
            per_species[species] = [identifier for identifier in members if identifier not in old_new_ids]

    if mode == "all_retained":
        additions = [row for row in audit_rows if row["retained_in_revised_bZIP_set"] == "Yes"]
    elif mode == "annotation_corrections_only":
        additions = [
            row
            for row in audit_rows
            if row["retained_in_revised_bZIP_set"] == "Yes"
            and row["independent_new_gene"] == "No"
        ]
    elif mode == "annotated_only":
        additions = []
    else:
        raise ValueError(mode)

    species_names = list(next(iter(groups.values())))
    for row in additions:
        orthogroup = row["orthogroup"]
        species_safe = row["species"].replace(".", "_")
        if orthogroup not in groups:
            raise RuntimeError(f"Audit references absent old orthogroup: {orthogroup}")
        if species_safe not in species_names:
            raise RuntimeError(f"Audit species absent from OrthoFinder columns: {species_safe}")
        groups[orthogroup][species_safe].append(row["deliverable_model_id"])

    return {orthogroup: group for orthogroup, group in groups.items() if group_members(group)}


def monte_carlo_fixed_margin_p(table: np.ndarray) -> tuple[float, int]:
    distribution = stats.random_table(
        table.sum(axis=1), table.sum(axis=0), seed=np.random.default_rng(SEED)
    )
    observed_log_probability = float(distribution.logpmf(table))
    extreme = 0
    completed = 0
    batch_size = 20_000
    while completed < SIMULATIONS:
        current = min(batch_size, SIMULATIONS - completed)
        sampled = distribution.rvs(size=current)
        log_probabilities = distribution.logpmf(sampled)
        extreme += int(np.count_nonzero(log_probabilities <= observed_log_probability + 1e-12))
        completed += current
    return float((extreme + 1) / (SIMULATIONS + 1)), extreme


def cramers_v(table: np.ndarray) -> tuple[float, float, int, float, float]:
    chi2, p_value, dof, expected = stats.chi2_contingency(table, correction=False)
    v = math.sqrt(chi2 / (table.sum() * min(table.shape[0] - 1, table.shape[1] - 1)))
    return float(chi2), float(p_value), int(dof), float(v), float(expected.min())


def duplication_analysis(
    label: str,
    groups: dict[str, dict[str, list[str]]],
    stats_table: pd.DataFrame,
    mapping_by_id: dict[str, dict[str, str]],
    duplication_lookup: dict[tuple[str, str], str],
    outdir: Path,
) -> tuple[dict[str, object], pd.DataFrame]:
    category_by_group = dict(zip(stats_table["Orthogroup"], stats_table["Category"]))
    assignment_rows = []
    missing_mapping = []
    missing_duplication = []
    for orthogroup, per_species in groups.items():
        for members in per_species.values():
            for identifier in members:
                mapping = mapping_by_id.get(identifier)
                if mapping is None:
                    missing_mapping.append(identifier)
                    continue
                key = (mapping["Species"], mapping["Original_ID"])
                duplication_type = duplication_lookup.get(key)
                if duplication_type is None:
                    missing_duplication.append(key)
                    continue
                assignment_rows.append(
                    {
                        "dataset": label,
                        "Orthogroup": orthogroup,
                        "Category": category_by_group[orthogroup],
                        "Species": mapping["Species"],
                        "OrthoFinder_ID": identifier,
                        "Original_ID": mapping["Original_ID"],
                        "Duplication_type": duplication_type,
                    }
                )
    if missing_mapping or missing_duplication:
        raise RuntimeError(
            f"Duplication join failed for {label}: missing mapping={missing_mapping[:5]}, "
            f"missing duplication={missing_duplication[:5]}"
        )

    assignments = pd.DataFrame(assignment_rows)
    row_order = [category for category in CATEGORY_ORDER if category in set(assignments["Category"])]
    full = (
        pd.crosstab(assignments["Category"], assignments["Duplication_type"])
        .reindex(index=row_order, columns=DUPLICATION_ORDER, fill_value=0)
        .astype(int)
    )
    full.to_csv(outdir / f"{label}_duplication_contingency_full.tsv", sep="\t")

    chi2, pearson_p, dof, v, expected_min = cramers_v(full.to_numpy())
    exact_p, extreme = monte_carlo_fixed_margin_p(full.to_numpy())
    collapsed = full.copy()
    collapsed["local/other"] = collapsed[["proximal", "tandem", "singleton"]].sum(axis=1)
    collapsed = collapsed[["WGD/segmental", "dispersed", "local/other"]]
    c_chi2, c_p, c_dof, c_v, c_expected_min = cramers_v(collapsed.to_numpy())
    collapsed.to_csv(outdir / f"{label}_duplication_contingency_collapsed.tsv", sep="\t")
    assignments.to_csv(outdir / f"{label}_duplication_gene_assignments.tsv", sep="\t", index=False)

    report: dict[str, object] = {
        "dataset": label,
        "genes_in_test": int(len(assignments)),
        "row_categories": row_order,
        "full_table": full.to_dict(orient="index"),
        "full_exact_test": {
            "test": "Fisher-Freeman-Halton Monte Carlo exact test",
            "simulations": SIMULATIONS,
            "seed": SEED,
            "p_value": exact_p,
            "extreme_simulations": extreme,
            "pearson_chi_square_for_effect_size": chi2,
            "df": dof,
            "cramers_v": v,
            "minimum_expected_count": expected_min,
            "pearson_p_value_not_primary": pearson_p,
        },
        "collapsed_test": {
            "test": "Pearson chi-square test",
            "chi_square": c_chi2,
            "df": c_dof,
            "p_value": c_p,
            "cramers_v": c_v,
            "minimum_expected_count": c_expected_min,
        },
    }
    return report, assignments


def main() -> None:
    args = parse_args()
    outdir = args.revision_root / "analysis"
    outdir.mkdir(exist_ok=True)
    source_code = args.revision_root / "source_code"
    source_code.mkdir(exist_ok=True)
    shutil.copy2(Path(__file__), source_code / Path(__file__).name)

    audit_rows = read_tsv(args.audit_tsv)
    old_mapping_rows = read_tsv(
        args.old_root / "tables" / "bZIP_orthofinder_id_mapping.tsv"
    )
    mapping_by_id = {row["OrthoFinder_ID"]: row for row in old_mapping_rows}
    old_new_ids = {
        row["OrthoFinder_ID"] for row in old_mapping_rows if row["Type"] == "new_predicted"
    }
    if len(old_new_ids) != 52:
        raise RuntimeError(f"Expected 52 old candidates, found {len(old_new_ids)}")

    old_results = args.old_root / "OrthoFinder_run" / "Results_Jun30"
    species, old_groups, old_unassigned = load_orthogroups(old_results)
    old_stats = group_stats(old_groups, species)

    old_metadata_path = args.duplication_tsv.parent / "ogg_copy_number_classification.tsv"
    old_metadata = pd.read_csv(old_metadata_path, sep="\t", dtype=str)
    old_meta_by_og = old_metadata.set_index("Orthogroup").to_dict(orient="index")
    old_stats["pangeneID"] = old_stats["Orthogroup"].map(
        lambda value: old_meta_by_og[value]["pangeneID"]
    )
    old_stats["Old_subfamily"] = old_stats["Orthogroup"].map(
        lambda value: old_meta_by_og[value]["Subfamily"]
    )
    old_stats.to_csv(outdir / "original_ogg_classification_recomputed.tsv", sep="\t", index=False)

    old_member_to_group = {
        identifier: orthogroup
        for orthogroup, group in old_groups.items()
        for identifier in group_members(group)
    }
    if len(old_member_to_group) != 1790:
        raise RuntimeError(f"Expected 1790 old assigned genes, found {len(old_member_to_group)}")

    dataset_modes = {
        "revised_audited_25": "all_retained",
        "sensitivity_no_independent_new": "annotation_corrections_only",
        "sensitivity_annotated_only": "annotated_only",
    }
    input_totals = {
        "revised_audited_25": 1764,
        "sensitivity_no_independent_new": 1743,
        "sensitivity_annotated_only": 1739,
    }

    summary_rows = [
        summary_row("original_1791", "original", old_stats, 1791, len(old_unassigned))
    ]
    fixed_groups: dict[str, dict[str, dict[str, list[str]]]] = {}
    fixed_stats: dict[str, pd.DataFrame] = {}
    de_novo_groups: dict[str, dict[str, dict[str, list[str]]]] = {}
    de_novo_stats: dict[str, pd.DataFrame] = {}

    for dataset, mode in dataset_modes.items():
        groups = fixed_membership_dataset(old_groups, old_new_ids, audit_rows, mode)
        table = group_stats(groups, species)
        table["pangeneID"] = table["Orthogroup"].map(
            lambda value: old_meta_by_og[value]["pangeneID"]
        )
        table.to_csv(outdir / f"fixed_membership_{dataset}_classification.tsv", sep="\t", index=False)
        fixed_groups[dataset] = groups
        fixed_stats[dataset] = table
        assigned = int(table["Total_genes"].sum())
        summary_rows.append(
            summary_row(dataset, "fixed_membership", table, input_totals[dataset], input_totals[dataset] - assigned)
        )

        results = result_dir(args.revision_root / "runs" / dataset)
        new_species, new_groups, new_unassigned = load_orthogroups(results)
        if new_species != species:
            raise RuntimeError(f"Species order changed for {dataset}")
        new_table = group_stats(new_groups, species)

        contributor_rows = []
        retained_old_group = {
            row["deliverable_model_id"]: row["orthogroup"]
            for row in audit_rows
            if row["retained_in_revised_bZIP_set"] == "Yes"
        }
        for orthogroup, group in new_groups.items():
            shared = Counter(
                old_member_to_group[identifier]
                for identifier in group_members(group)
                if identifier in old_member_to_group
            )
            added = Counter(
                retained_old_group[identifier]
                for identifier in group_members(group)
                if identifier in retained_old_group
            )
            combined = sorted(set(shared) | set(added))
            contributor_rows.append(
                {
                    "dataset": dataset,
                    "New_Orthogroup": orthogroup,
                    "Old_orthogroup_contributors": ";".join(combined) or "NA",
                    "Old_pangene_contributors": ";".join(
                        old_meta_by_og[value]["pangeneID"] for value in combined
                    ) or "NA",
                    "Shared_old_member_counts": ";".join(
                        f"{value}:{shared[value]}" for value in sorted(shared)
                    ) or "NA",
                    "Added_final_model_counts": ";".join(
                        f"{value}:{added[value]}" for value in sorted(added)
                    ) or "NA",
                    "Merge_of_multiple_old_OGGs": "Yes" if len(combined) > 1 else "No",
                }
            )
        contributors = pd.DataFrame(contributor_rows)
        new_table = new_table.merge(
            contributors, left_on="Orthogroup", right_on="New_Orthogroup", how="left"
        ).drop(columns=["New_Orthogroup"])
        new_table.to_csv(outdir / f"de_novo_{dataset}_classification.tsv", sep="\t", index=False)
        contributors.to_csv(outdir / f"de_novo_{dataset}_old_group_contributors.tsv", sep="\t", index=False)
        de_novo_groups[dataset] = new_groups
        de_novo_stats[dataset] = new_table
        summary_rows.append(
            summary_row(dataset, "de_novo_orthofinder", new_table, input_totals[dataset], len(new_unassigned))
        )

    summary = pd.DataFrame(summary_rows)
    summary.to_csv(outdir / "ogg_revision_summary.tsv", sep="\t", index=False)

    mapping_rows = []
    change_rows = []
    old_category = dict(zip(old_stats["Orthogroup"], old_stats["Category"]))
    for dataset, groups in de_novo_groups.items():
        new_member_to_group = {
            identifier: orthogroup
            for orthogroup, group in groups.items()
            for identifier in group_members(group)
        }
        new_category = dict(
            zip(de_novo_stats[dataset]["Orthogroup"], de_novo_stats[dataset]["Category"])
        )
        for old_group, group in old_groups.items():
            surviving_annotated = sorted(
                identifier
                for identifier in group_members(group)
                if identifier not in old_new_ids
            )
            targets = Counter(
                new_member_to_group[identifier]
                for identifier in surviving_annotated
                if identifier in new_member_to_group
            )
            target_names = sorted(targets)
            mapping_rows.append(
                {
                    "dataset": dataset,
                    "Old_Orthogroup": old_group,
                    "Old_pangeneID": old_meta_by_og[old_group]["pangeneID"],
                    "Old_category": old_category[old_group],
                    "Surviving_annotated_members": len(surviving_annotated),
                    "Mapped_new_OGGs": ";".join(target_names) or "NA",
                    "Mapped_member_counts": ";".join(
                        f"{value}:{targets[value]}" for value in target_names
                    ) or "NA",
                    "Mapping_status": (
                        "one_to_one" if len(target_names) == 1 else "split" if len(target_names) > 1 else "lost"
                    ),
                }
            )
            if len(target_names) == 1:
                current = new_category[target_names[0]]
                change_rows.append(
                    {
                        "dataset": dataset,
                        "Old_Orthogroup": old_group,
                        "Old_pangeneID": old_meta_by_og[old_group]["pangeneID"],
                        "Old_category": old_category[old_group],
                        "New_Orthogroup": target_names[0],
                        "New_category": current,
                        "Category_changed": "Yes" if current != old_category[old_group] else "No",
                    }
                )
    pd.DataFrame(mapping_rows).to_csv(outdir / "old_to_de_novo_ogg_mapping.tsv", sep="\t", index=False)
    changes = pd.DataFrame(change_rows)
    changes.to_csv(outdir / "old_to_de_novo_category_changes.tsv", sep="\t", index=False)

    fixed_change_rows = []
    for dataset, table in fixed_stats.items():
        current_category = dict(zip(table["Orthogroup"], table["Category"]))
        for old_group in old_groups:
            new_value = current_category.get(old_group, "absent")
            fixed_change_rows.append(
                {
                    "dataset": dataset,
                    "Old_Orthogroup": old_group,
                    "Old_pangeneID": old_meta_by_og[old_group]["pangeneID"],
                    "Old_category": old_category[old_group],
                    "New_category": new_value,
                    "Category_changed": "Yes" if new_value != old_category[old_group] else "No",
                }
            )
    fixed_changes = pd.DataFrame(fixed_change_rows)
    fixed_changes.to_csv(outdir / "fixed_membership_category_changes.tsv", sep="\t", index=False)

    duplication_rows = read_tsv(args.duplication_tsv)
    duplication_lookup = {
        (row["species"], row["geneID"]): row["type"] for row in duplication_rows
    }
    if len(duplication_lookup) != len(duplication_rows):
        raise RuntimeError("Duplication lookup contains duplicate species/geneID keys")

    fixed_dup_report, _ = duplication_analysis(
        "fixed_annotated_only",
        fixed_groups["sensitivity_annotated_only"],
        fixed_stats["sensitivity_annotated_only"],
        mapping_by_id,
        duplication_lookup,
        outdir,
    )
    denovo_dup_report, _ = duplication_analysis(
        "de_novo_annotated_only",
        de_novo_groups["sensitivity_annotated_only"],
        de_novo_stats["sensitivity_annotated_only"],
        mapping_by_id,
        duplication_lookup,
        outdir,
    )

    cr021_groups = [
        orthogroup
        for orthogroup, meta in old_meta_by_og.items()
        if meta["pangeneID"] == "PtbZIP.CR021"
    ]
    if len(cr021_groups) != 1:
        raise RuntimeError(f"Expected one PtbZIP.CR021 group, found {cr021_groups}")
    cr021 = cr021_groups[0]
    cr021_rows = []
    for _, row in summary.iterrows():
        dataset = row["dataset"]
        mode = row["analysis_mode"]
        if mode == "original":
            table = old_stats
            groups = old_groups
            target = cr021
        elif mode == "fixed_membership":
            table = fixed_stats[dataset]
            groups = fixed_groups[dataset]
            target = cr021
        else:
            mapping_subset = pd.DataFrame(mapping_rows)
            mapping_subset = mapping_subset[
                (mapping_subset["dataset"] == dataset)
                & (mapping_subset["Old_Orthogroup"] == cr021)
            ]
            mapped = mapping_subset.iloc[0]["Mapped_new_OGGs"]
            if mapped == "NA" or ";" in mapped:
                continue
            target = mapped
            table = de_novo_stats[dataset]
            groups = de_novo_groups[dataset]
        stat_row = table[table["Orthogroup"] == target].iloc[0]
        cr021_rows.append(
            {
                "dataset": dataset,
                "analysis_mode": mode,
                "Orthogroup": target,
                "Category": stat_row["Category"],
                "Species_present": int(stat_row["Species_present"]),
                "Total_genes": int(stat_row["Total_genes"]),
                "Copy_number_variable": bool(stat_row["Copy_number_variable"]),
            }
        )
    pd.DataFrame(cr021_rows).to_csv(outdir / "PtbZIP_CR021_sensitivity.tsv", sep="\t", index=False)

    fixed_strict_summary = summary[
        (summary["dataset"] == "sensitivity_annotated_only")
        & (summary["analysis_mode"] == "fixed_membership")
    ].iloc[0].to_dict()
    denovo_strict_summary = summary[
        (summary["dataset"] == "sensitivity_annotated_only")
        & (summary["analysis_mode"] == "de_novo_orthofinder")
    ].iloc[0].to_dict()
    main_summary = summary[
        (summary["dataset"] == "revised_audited_25")
        & (summary["analysis_mode"] == "de_novo_orthofinder")
    ].iloc[0].to_dict()

    report = {
        "control_reproduced_original": True,
        "original": summary.iloc[0].to_dict(),
        "fixed_membership_strict_exclusion": fixed_strict_summary,
        "de_novo_strict_exclusion": denovo_strict_summary,
        "de_novo_revised_main": main_summary,
        "fixed_membership_changed_pangenes": fixed_changes[
            (fixed_changes["dataset"] == "sensitivity_annotated_only")
            & (fixed_changes["Category_changed"] == "Yes")
        ][["Old_Orthogroup", "Old_pangeneID", "Old_category", "New_category"]].to_dict(
            orient="records"
        ),
        "de_novo_strict_changed_pangenes": changes[
            (changes["dataset"] == "sensitivity_annotated_only")
            & (changes["Category_changed"] == "Yes")
        ][["Old_Orthogroup", "Old_pangeneID", "Old_category", "New_Orthogroup", "New_category"]].to_dict(
            orient="records"
        ),
        "PtbZIP_CR021": cr021_rows,
        "duplication_tests": {
            "fixed_membership_strict_exclusion": fixed_dup_report,
            "de_novo_strict_exclusion": denovo_dup_report,
        },
    }
    (outdir / "reviewer_key_findings.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2, default=str))


if __name__ == "__main__":
    main()
