#!/usr/bin/env python3
"""Re-map validated bZIP genes to revised OGGs and summarize duplication modes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
from scipy import stats


CLASS_ORDER = ["core", "softcore", "shell", "cloud"]
MODE_ORDER = ["WGD/segmental", "dispersed", "proximal", "tandem", "singleton"]
SIMULATIONS = 1_000_000
SEED = 20260721

# Annotation-corrected models replace existing gene models and therefore inherit
# the genome-wide MCScanX class of the corrected locus, not the provisional
# BITACORA fragment that was appended during the historical analysis.
ANNOTATION_CORRECTION_MODES = {
    ("Populus_alba", "CM100082.1_exon1_exon4_exon6"): "tandem",
    ("Populus_euphratica", "GWHAAYU00000007_exon13"): "tandem",
    ("Populus_pruinosa", "Chr7_exon13_exon12_exon25_exon1_exon7_exon11"): "tandem",
    ("Populus_tremula", "chr9_exon1"): "singleton",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s4-json", required=True, type=Path)
    parser.add_argument("--old-s7-json", required=True, type=Path)
    parser.add_argument("--candidate-audit-json", required=True, type=Path)
    parser.add_argument("--annotated-duplication-tsv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def normalize_species(value: str) -> str:
    return "Populus_deltoides" if value == "Populus_deltoide" else value


def normalize_class(value: str) -> str:
    return value.replace("soft-core", "softcore").replace("private", "cloud")


def parse_s4(path: Path) -> list[dict[str, str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload["values"]
    header = rows[1]
    index = {str(value): i for i, value in enumerate(header) if value not in (None, "")}
    members: list[dict[str, str]] = []
    for row in rows[2:81]:
        if not row or row[0] in (None, ""):
            continue
        orthogroup = str(row[index["Orthogroup"]])
        pangene = str(row[index["Pangene_ID"]])
        category = normalize_class(str(row[index["Classification"]]))
        cell = str(row[index["All_genes_in_OGG"]])
        for entry in cell.split(", "):
            species, gene_id, source = entry.split("|", 2)
            members.append(
                {
                    "orthogroup": orthogroup,
                    "class": category,
                    "pangeneID": pangene,
                    "geneID": gene_id,
                    "source": source,
                    "species": normalize_species(species),
                }
            )
    return members


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_old_s7(path: Path) -> list[dict[str, str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload["values"]
    records = []
    for row in rows[2:1792]:
        if len(row) < 5 or row[2] in (None, ""):
            continue
        records.append(
            {
                "class": normalize_class(str(row[0])),
                "pangeneID": str(row[1]),
                "geneID": str(row[2]),
                "type": str(row[3]),
                "species": normalize_species(str(row[4])),
            }
        )
    return records


def parse_audit(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload["values"]
    header = [str(value) for value in rows[0]]
    result = {}
    for values in rows[1:]:
        row = {header[i]: values[i] if i < len(values) else None for i in range(len(header))}
        if row["Retained in revised bZIP set"] != "Yes":
            continue
        key = (normalize_species(str(row["Species"])), str(row["Candidate ID"]))
        result[key] = row
    return result


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


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
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    members = parse_s4(args.s4_json)
    old_s7 = parse_old_s7(args.old_s7_json)
    audit = parse_audit(args.candidate_audit_json)
    annotated = read_tsv(args.annotated_duplication_tsv)

    annotated_lookup = {
        (normalize_species(row["species"]), row["geneID"]): row["type"] for row in annotated
    }
    candidate_lookup = {
        (row["species"], row["geneID"]): row["type"] for row in old_s7
    }

    assignments: list[dict[str, object]] = []
    missing: list[dict[str, str]] = []
    for member in members:
        key = (member["species"], member["geneID"])
        if member["source"] == "annotated":
            mode = annotated_lookup.get(key)
            source = "existing genome-wide MCScanX call for unchanged annotated model"
        elif member["source"] == "annotation_corrected":
            mode = ANNOTATION_CORRECTION_MODES.get(key)
            source = "genome-wide MCScanX call inherited from corrected annotated locus"
        else:
            mode = candidate_lookup.get(key)
            source = "existing genome-wide MCScanX call for audited candidate locus"
        if mode is None:
            missing.append(member)
            continue
        assignments.append(
            {
                **member,
                "type": mode,
                "duplication_call_source": source,
            }
        )

    if missing:
        raise RuntimeError(f"Missing duplication calls for {len(missing)} members: {missing[:8]}")
    if len(assignments) != 1762:
        raise RuntimeError(f"Expected 1,762 assigned proteins, found {len(assignments)}")
    if len({(row['species'], row['geneID']) for row in assignments}) != 1762:
        raise RuntimeError("Current gene keys are not unique")

    retained_rows = []
    for row in assignments:
        if row["source"] == "annotated":
            continue
        audit_row = audit[(row["species"], row["geneID"])]
        retained_rows.append(
            {
                "species": row["species"],
                "geneID": row["geneID"],
                "final_classification": audit_row["Final classification"],
                "original_sequence": audit_row["Original genomic sequence"],
                "original_start": audit_row["Original model start"],
                "original_end": audit_row["Original model end"],
                "selected_sequence": audit_row["Selected/recovered genomic sequence"],
                "selected_start": audit_row["Selected/recovered model start"],
                "selected_end": audit_row["Selected/recovered model end"],
                "same_genomic_sequence": "Yes"
                if audit_row["Original genomic sequence"] == audit_row["Selected/recovered genomic sequence"]
                else "No",
                "duplication_type": row["type"],
                "orthogroup": row["orthogroup"],
                "pangeneID": row["pangeneID"],
            }
        )

    write_tsv(
        args.output_dir / "revised_duplication_gene_assignments.tsv",
        assignments,
        [
            "class",
            "pangeneID",
            "geneID",
            "type",
            "species",
            "orthogroup",
            "source",
            "duplication_call_source",
        ],
    )
    write_tsv(
        args.output_dir / "retained_model_duplication_audit.tsv",
        retained_rows,
        list(retained_rows[0]),
    )

    full = np.array(
        [
            [
                sum(1 for row in assignments if row["class"] == category and row["type"] == mode)
                for mode in MODE_ORDER
            ]
            for category in CLASS_ORDER
        ],
        dtype=int,
    )
    exact, extreme = exact_p(full)
    primary = {
        "test": "Fisher-Freeman-Halton Monte Carlo exact test",
        "simulations": SIMULATIONS,
        "seed": SEED,
        "p_value": exact,
        "extreme_simulations": extreme,
        **metrics(full),
    }

    comparable_rows = []
    for row in assignments:
        comparable_rows.append(
            {
                **row,
                "comparable_class": "shell" if row["class"] == "cloud" else row["class"],
                "comparable_mode": row["type"]
                if row["type"] in {"WGD/segmental", "dispersed"}
                else "local/other",
            }
        )
    comparable_classes = ["core", "softcore", "shell"]
    comparable_modes = ["WGD/segmental", "dispersed", "local/other"]
    collapsed = np.array(
        [
            [
                sum(
                    1
                    for row in comparable_rows
                    if row["comparable_class"] == category
                    and row["comparable_mode"] == mode
                )
                for mode in comparable_modes
            ]
            for category in comparable_classes
        ],
        dtype=int,
    )
    report = {
        "genes": len(assignments),
        "source_counts": Counter(row["source"] for row in assignments),
        "mode_counts": Counter(row["type"] for row in assignments),
        "class_counts": Counter(row["class"] for row in assignments),
        "full_table": {
            category: {mode: int(full[i, j]) for j, mode in enumerate(MODE_ORDER)}
            for i, category in enumerate(CLASS_ORDER)
        },
        "full_test": primary,
        "collapsed_table": {
            category: {mode: int(collapsed[i, j]) for j, mode in enumerate(comparable_modes)}
            for i, category in enumerate(comparable_classes)
        },
        "collapsed_test": {"test": "Pearson chi-square test", **metrics(collapsed)},
    }
    (args.output_dir / "revised_duplication_statistics.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, default=dict) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, ensure_ascii=False, indent=2, default=dict))


if __name__ == "__main__":
    main()
