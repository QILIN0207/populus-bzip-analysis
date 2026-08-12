#!/usr/bin/env python3
"""Matched non-bZIP TE background and OGG-level effect-size analysis.

One coding transcript is selected per gene. Controls are matched within species
and chromosome by transcript length and relative sequence position.
"""

from __future__ import annotations

import argparse
import bisect
import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from te_io import (
    EXCLUDED_REPEAT_CLASSES,
    TERecord,
    merge_coverage,
    read_fai,
    read_repeatmasker,
    read_seq_mapping,
)

ANCHOR_RE = re.compile(r"^Chr(?:[1-9]|1[0-9])$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bzip-summary", type=Path, required=True)
    parser.add_argument("--ogg-metrics", type=Path, required=True)
    parser.add_argument("--te-root", type=Path, required=True)
    parser.add_argument("--mcscanx-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--matching-replicates", type=int, default=1000)
    parser.add_argument("--bootstrap-replicates", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260810)
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.strip().strip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def internal_species(display: str) -> str:
    if display == "Populus_deltoides":
        return "Populus_deltoide"
    if display in {"Populus_alba_var_pyramidalis", "Populus_alba_var__pyramidalis"}:
        return "Populus_alba_var._pyramidalis"
    return display


def display_species(internal: str) -> str:
    return "Populus_deltoides" if internal == "Populus_deltoide" else internal


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    preferred: bool


def read_coding_transcripts(gff: Path, seq_map: dict[str, str], fai: dict[str, int]) -> pd.DataFrame:
    transcripts: dict[str, Transcript] = {}
    coding: set[str] = set()
    for line in gff.open(encoding="utf-8", errors="replace"):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        feature = fields[2]
        attrs = parse_attributes(fields[8])
        if feature in {"mRNA", "transcript"}:
            transcript_id = attrs.get("ID", "").split(",", 1)[0]
            if not transcript_id:
                continue
            parents = attrs.get("Parent", "").split(",")
            gene_id = parents[0] if parents and parents[0] else transcript_id
            chrom = seq_map.get(fields[0], fields[0])
            if chrom not in fai:
                continue
            preferred = any(
                attrs.get(key, "").lower() in {"1", "true", "yes"}
                for key in ("longest", "representative", "primary")
            )
            transcripts[transcript_id] = Transcript(
                transcript_id, gene_id, chrom, int(fields[3]), int(fields[4]), fields[6], preferred
            )
        elif feature == "CDS":
            coding.update(parent for parent in attrs.get("Parent", "").split(",") if parent)

    by_gene: dict[str, list[Transcript]] = defaultdict(list)
    for transcript_id, record in transcripts.items():
        if transcript_id in coding:
            by_gene[record.gene_id].append(record)
    rows: list[dict[str, object]] = []
    for gene_id, records in by_gene.items():
        records.sort(key=lambda r: (-int(r.preferred), -(r.end - r.start + 1), r.transcript_id))
        r = records[0]
        seq_len = fai[r.chrom]
        rows.append(
            {
                "control_gene_id": gene_id,
                "control_transcript_id": r.transcript_id,
                "gene_chr": r.chrom,
                "gene_start": r.start,
                "gene_end": r.end,
                "gene_strand": r.strand,
                "gene_length_bp": r.end - r.start + 1,
                "sequence_length_bp": seq_len,
                "relative_midpoint": (((r.start - 1) + r.end) / 2.0) / seq_len,
                "anchored_chromosome": bool(ANCHOR_RE.match(r.chrom)),
            }
        )
    return pd.DataFrame(rows)


def remove_bzip_loci(controls: pd.DataFrame, bzip: pd.DataFrame) -> pd.DataFrame:
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for rec in bzip.to_dict("records"):
        by_chrom[str(rec["gene_chr"])].append((int(rec["gene_start"]), int(rec["gene_end"])))
    keep: list[bool] = []
    for rec in controls.to_dict("records"):
        start, end = int(rec["gene_start"]), int(rec["gene_end"])
        keep.append(not any(max(start, a) <= min(end, b) for a, b in by_chrom.get(str(rec["gene_chr"]), [])))
    return controls.loc[keep].reset_index(drop=True)


def calculate_window_te(frame: pd.DataFrame, fai: dict[str, int], tes: dict[str, list[TERecord]]) -> pd.DataFrame:
    starts = {chrom: [r.start0 for r in records] for chrom, records in tes.items()}
    covered: list[int] = []
    lengths: list[int] = []
    for rec in frame.to_dict("records"):
        chrom = str(rec["gene_chr"])
        start0, end0 = int(rec["gene_start"]) - 1, int(rec["gene_end"])
        window_start, window_end = max(0, start0 - 2000), min(fai[chrom], end0 + 2000)
        intervals: list[tuple[int, int]] = []
        values = tes.get(chrom, [])
        stop = bisect.bisect_left(starts.get(chrom, []), window_end)
        for te in values[:stop]:
            if te.end0 <= window_start:
                continue
            left, right = max(window_start, te.start0), min(window_end, te.end0)
            if right > left:
                intervals.append((left, right))
        covered.append(merge_coverage(intervals))
        lengths.append(window_end - window_start)
    result = frame.copy()
    result["te_covered_bp"] = covered
    result["window_bp"] = lengths
    result["te_coverage_fraction"] = result["te_covered_bp"] / result["window_bp"]
    result["TE_positive"] = result["te_covered_bp"] > 0
    return result


def read_te_library_summary(path: Path) -> tuple[int, int]:
    total = retained = 0
    for line in path.open(encoding="utf-8", errors="replace"):
        if line.startswith(">"):
            total += 1
            header = line[1:].strip().split()[0]
            repeat_class = header.split("#", 1)[1] if "#" in header else "Unknown"
            retained += int(repeat_class not in EXCLUDED_REPEAT_CLASSES)
    return total, retained


def genome_te_summary(species: str, fai: dict[str, int], tes: dict[str, list[TERecord]], library: Path) -> dict[str, object]:
    assembly_bp = sum(fai.values())
    covered_bp = sum(merge_coverage((r.start0, r.end0) for r in values) for values in tes.values())
    total_families, retained_families = read_te_library_summary(library)
    return {
        "Species": display_species(species),
        "assembly_sequences_n": len(fai),
        "assembly_bp": assembly_bp,
        "filtered_TE_hits_n": sum(len(values) for values in tes.values()),
        "nonredundant_TE_covered_bp": covered_bp,
        "TE_assembly_fraction": covered_bp / assembly_bp,
        "HiTE_consensus_families_all_n": total_families,
        "HiTE_consensus_families_analyzed_n": retained_families,
        "repeat_filter": "Excluded Simple_repeat, Low_complexity, Satellite, rRNA, scRNA, snRNA, srpRNA, tRNA, and RNA",
    }


def candidate_distance(target: pd.Series, controls: pd.DataFrame) -> np.ndarray:
    length_term = np.abs(np.log2(controls["gene_length_bp"].to_numpy(float) / float(target["gene_length_bp"]))) / 0.25
    position_term = np.abs(controls["relative_midpoint"].to_numpy(float) - float(target["relative_midpoint"])) / 0.05
    distance = length_term + position_term
    if not bool(target["anchored_chromosome"]):
        distance += np.abs(np.log2(controls["sequence_length_bp"].to_numpy(float) / float(target["sequence_length_bp"]))) / 0.5
    return distance


def build_match_lists(bzip: pd.DataFrame, controls: pd.DataFrame) -> list[dict[str, object]]:
    anchored = controls[controls["anchored_chromosome"]]
    unplaced = controls[~controls["anchored_chromosome"]]
    result: list[dict[str, object]] = []
    for target_index, target in bzip.reset_index(drop=True).iterrows():
        if bool(target["anchored_chromosome"]):
            pool, scope = anchored[anchored["gene_chr"] == target["gene_chr"]], "same chromosome"
        else:
            same = unplaced[unplaced["gene_chr"] == target["gene_chr"]]
            pool, scope = (same, "same unplaced scaffold") if len(same) >= 5 else (unplaced, "unplaced-scaffold pool")
        if pool.empty:
            raise ValueError(f"No matched control pool for {target['OrthoFinder_ID']} on {target['gene_chr']}")
        distances = candidate_distance(target, pool)
        order = np.argsort(distances, kind="mergesort")
        result.append(
            {
                "target_index": target_index,
                "target_id": target["OrthoFinder_ID"],
                "scope": scope,
                "control_indices": pool.index.to_numpy(int)[order],
                "distances": distances[order],
            }
        )
    return result


def empirical_two_sided(reference: float, samples: np.ndarray) -> float:
    upper = (np.count_nonzero(samples >= reference) + 1) / (len(samples) + 1)
    lower = (np.count_nonzero(samples <= reference) + 1) / (len(samples) + 1)
    return min(1.0, 2.0 * min(upper, lower))


def matched_background(
    bzip: pd.DataFrame, controls: pd.DataFrame, replicates: int, seed: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    bzip, controls = bzip.reset_index(drop=True).copy(), controls.reset_index(drop=True).copy()
    match_lists = build_match_lists(bzip, controls)
    species = sorted(bzip["Species"].unique())
    species_index = {s: i for i, s in enumerate(species)}
    rep_rows: list[dict[str, object]] = []
    per_species = np.zeros((replicates, len(species)), dtype=int)
    first_mapping: list[dict[str, object]] = []
    for replicate in range(replicates):
        used: set[int] = set()
        selected: dict[int, tuple[int, float, str]] = {}
        for target_index in rng.permutation(len(match_lists)):
            entry = match_lists[int(target_index)]
            available = [i for i, control_index in enumerate(entry["control_indices"]) if int(control_index) not in used]
            if not available:
                raise ValueError(f"Control pool exhausted for {entry['target_id']}")
            chosen_position = int(rng.choice(available[: min(5, len(available))]))
            control_index = int(entry["control_indices"][chosen_position])
            used.add(control_index)
            selected[int(entry["target_index"])] = (
                control_index, float(entry["distances"][chosen_position]), str(entry["scope"])
            )
        ordered = [selected[i] for i in range(len(bzip))]
        chosen = controls.loc[[x[0] for x in ordered]].reset_index(drop=True)
        te_positive = chosen["TE_positive"].to_numpy(bool)
        for i, sp in enumerate(bzip["Species"]):
            per_species[replicate, species_index[str(sp)]] += int(te_positive[i])
        log_ratio = np.abs(np.log2(chosen["gene_length_bp"].to_numpy(float) / bzip["gene_length_bp"].to_numpy(float)))
        pos_diff = np.abs(chosen["relative_midpoint"].to_numpy(float) - bzip["relative_midpoint"].to_numpy(float))
        rep_rows.append(
            {
                "replicate": replicate + 1,
                "matched_control_genes_n": len(chosen),
                "TE_positive_control_genes_n": int(te_positive.sum()),
                "TE_positive_control_fraction": float(te_positive.mean()),
                "aggregate_TE_coverage_fraction": float(chosen["te_covered_bp"].sum() / chosen["window_bp"].sum()),
                "median_absolute_log2_length_ratio": float(np.median(log_ratio)),
                "median_absolute_relative_position_difference": float(np.median(pos_diff)),
            }
        )
        if replicate == 0:
            for i, (control_index, distance, scope) in enumerate(ordered):
                target, control = bzip.iloc[i], controls.loc[control_index]
                first_mapping.append(
                    {
                        "Species": target["Species"],
                        "bZIP_OrthoFinder_ID": target["OrthoFinder_ID"],
                        "bZIP_chr": target["gene_chr"],
                        "bZIP_start": target["gene_start"],
                        "bZIP_end": target["gene_end"],
                        "bZIP_length_bp": target["gene_length_bp"],
                        "bZIP_relative_midpoint": target["relative_midpoint"],
                        "bZIP_TE_positive": target["TE_positive"],
                        "control_gene_id": control["control_gene_id"],
                        "control_transcript_id": control["control_transcript_id"],
                        "control_chr": control["gene_chr"],
                        "control_start": control["gene_start"],
                        "control_end": control["gene_end"],
                        "control_length_bp": control["gene_length_bp"],
                        "control_relative_midpoint": control["relative_midpoint"],
                        "control_TE_positive": control["TE_positive"],
                        "match_scope": scope,
                        "match_distance": distance,
                    }
                )

    replicate_frame = pd.DataFrame(rep_rows)
    bzip_positive, bzip_fraction = int(bzip["TE_positive"].sum()), float(bzip["TE_positive"].mean())
    control_fraction = replicate_frame["TE_positive_control_fraction"].to_numpy(float)
    difference = bzip_fraction - control_fraction
    overall = pd.DataFrame(
        [{
            "comparison": "All 19 genomes: bZIP vs matched non-bZIP genes",
            "bZIP_genes_n": len(bzip),
            "bZIP_TE_positive_n": bzip_positive,
            "bZIP_TE_positive_fraction": bzip_fraction,
            "matching_replicates": replicates,
            "matched_control_genes_per_replicate": len(bzip),
            "matched_control_TE_positive_fraction_mean": float(control_fraction.mean()),
            "matched_control_TE_positive_fraction_median": float(np.median(control_fraction)),
            "matched_control_TE_positive_fraction_ci_low": float(np.quantile(control_fraction, 0.025)),
            "matched_control_TE_positive_fraction_ci_high": float(np.quantile(control_fraction, 0.975)),
            "bZIP_minus_control_fraction_mean": float(difference.mean()),
            "bZIP_minus_control_fraction_ci_low": float(np.quantile(difference, 0.025)),
            "bZIP_minus_control_fraction_ci_high": float(np.quantile(difference, 0.975)),
            "empirical_two_sided_p": empirical_two_sided(bzip_fraction, control_fraction),
            "seed": seed,
        }]
    )
    by_species_rows: list[dict[str, object]] = []
    for sp in species:
        subset = bzip[bzip["Species"] == sp]
        n, bpos = len(subset), int(subset["TE_positive"].sum())
        cf = per_species[:, species_index[sp]] / n
        diff = bpos / n - cf
        by_species_rows.append(
            {
                "Species": sp,
                "bZIP_genes_n": n,
                "bZIP_TE_positive_n": bpos,
                "bZIP_TE_positive_fraction": bpos / n,
                "matched_control_TE_positive_fraction_mean": float(cf.mean()),
                "matched_control_TE_positive_fraction_ci_low": float(np.quantile(cf, 0.025)),
                "matched_control_TE_positive_fraction_ci_high": float(np.quantile(cf, 0.975)),
                "bZIP_minus_control_fraction_mean": float(diff.mean()),
                "bZIP_minus_control_fraction_ci_low": float(np.quantile(diff, 0.025)),
                "bZIP_minus_control_fraction_ci_high": float(np.quantile(diff, 0.975)),
                "empirical_two_sided_p": empirical_two_sided(bpos / n, cf),
            }
        )
    return overall, pd.DataFrame(by_species_rows), replicate_frame, pd.DataFrame(first_mapping)


def cliff_delta(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.sign(x[:, None] - y[None, :]).mean())


def cliff_bootstrap_ci(
    x: np.ndarray, y: np.ndarray, replicates: int, rng: np.random.Generator
) -> tuple[float, float, float]:
    observed = cliff_delta(x, y)
    values = np.empty(replicates, dtype=float)
    for start in range(0, replicates, 500):
        stop = min(start + 500, replicates)
        count = stop - start
        xb = x[rng.integers(0, len(x), size=(count, len(x)))]
        yb = y[rng.integers(0, len(y), size=(count, len(y)))]
        values[start:stop] = np.sign(xb[:, :, None] - yb[:, None, :]).mean(axis=(1, 2))
    return observed, float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975))


def build_cliff_table(ogg: pd.DataFrame, replicates: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    is_variable = ogg["Copy_number_variable"].astype(str).str.lower() == "true"
    variable, invariant = ogg[is_variable], ogg[~is_variable]
    metrics = [
        ("Aggregate TE coverage fraction", "Aggregate_TE_coverage_fraction"),
        ("TE-positive-gene fraction", "TE_positive_fraction"),
        ("Mean gene TE coverage fraction", "Mean_gene_TE_coverage_fraction"),
    ]
    rows: list[dict[str, object]] = []
    for label, column in metrics:
        x, y = variable[column].to_numpy(float), invariant[column].to_numpy(float)
        statistic, p_value = mannwhitneyu(x, y, alternative="two-sided")
        delta, low, high = cliff_bootstrap_ci(x, y, replicates, rng)
        rows.append(
            {
                "outcome": label,
                "unit": "OGG",
                "CNV_OGGs_n": len(x),
                "invariant_OGGs_n": len(y),
                "test": "Two-sided Mann-Whitney U",
                "U_statistic": statistic,
                "p_value": p_value,
                "effect_size": "Cliff's delta",
                "Cliffs_delta": delta,
                "bootstrap_CI_level": 0.95,
                "Cliffs_delta_CI_low": low,
                "Cliffs_delta_CI_high": high,
                "bootstrap_replicates": replicates,
                "bootstrap_unit": "OGG; groups resampled independently",
                "seed": seed,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=False)
    bzip = pd.read_csv(args.bzip_summary, sep="\t")
    ogg = pd.read_csv(args.ogg_metrics, sep="\t")
    required = {"Species", "OrthoFinder_ID", "gene_chr", "gene_start", "gene_end", "TE_positive"}
    if not required.issubset(bzip.columns):
        raise ValueError(f"Missing bZIP columns: {sorted(required - set(bzip.columns))}")
    if len(bzip) != 1764 or bzip["OrthoFinder_ID"].nunique() != 1764:
        raise ValueError("Expected 1,764 unique bZIP proteins")
    bzip["TE_positive"] = bzip["TE_positive"].astype(str).str.lower() == "true"
    bzip["gene_length_bp"] = bzip["gene_end"].astype(int) - bzip["gene_start"].astype(int) + 1

    control_frames: list[pd.DataFrame] = []
    genome_rows: list[dict[str, object]] = []
    pool_rows: list[dict[str, object]] = []
    bzip_frames: list[pd.DataFrame] = []
    for species_display in sorted(bzip["Species"].unique()):
        species = internal_species(str(species_display))
        te_dir = args.te_root / species
        gff = args.mcscanx_root / species / f"{species}.gene.original.gff"
        fai_path = te_dir / f"{species}.full_for_HiTE.fa.fai"
        repeat_path = te_dir / "S12_TE2kb_from_confTE_fullgenome" / f"{species}.full_for_HiTE.fa.out"
        library = te_dir / "HiTE_out_full_genome" / "confident_TE.cons.fa"
        for source in (gff, fai_path, repeat_path, library):
            if not source.exists():
                raise FileNotFoundError(source)
        fai = read_fai(fai_path)
        seq_map = read_seq_mapping(te_dir / f"{species}.full_chr_mapping.txt")
        controls = read_coding_transcripts(gff, seq_map, fai)
        species_bzip = bzip[bzip["Species"] == species_display].copy()
        species_bzip["sequence_length_bp"] = species_bzip["gene_chr"].map(fai)
        if species_bzip["sequence_length_bp"].isna().any():
            raise ValueError(f"Missing sequence length for bZIP genes in {species_display}")
        species_bzip["relative_midpoint"] = (
            ((species_bzip["gene_start"].astype(int) - 1) + species_bzip["gene_end"].astype(int)) / 2.0
        ) / species_bzip["sequence_length_bp"]
        species_bzip["anchored_chromosome"] = species_bzip["gene_chr"].astype(str).str.match(ANCHOR_RE)
        controls_before = len(controls)
        controls = remove_bzip_loci(controls, species_bzip)
        tes = read_repeatmasker(repeat_path, set(fai))
        controls = calculate_window_te(controls, fai, tes)
        controls.insert(0, "Species", species_display)
        control_frames.append(controls)
        bzip_frames.append(species_bzip)
        genome_rows.append(genome_te_summary(species, fai, tes, library))
        pool_rows.append(
            {
                "Species": species_display,
                "representative_coding_genes_before_bZIP_exclusion": controls_before,
                "eligible_non_bZIP_control_genes": len(controls),
                "control_pool_TE_positive_genes": int(controls["TE_positive"].sum()),
                "control_pool_TE_positive_fraction": float(controls["TE_positive"].mean()),
                "anchored_control_genes": int(controls["anchored_chromosome"].sum()),
                "unplaced_scaffold_control_genes": int((~controls["anchored_chromosome"]).sum()),
            }
        )

    all_bzip = pd.concat(bzip_frames, ignore_index=True)
    all_controls = pd.concat(control_frames, ignore_index=True)
    overall, by_species, replicates, first_mapping = matched_background(
        all_bzip, all_controls, args.matching_replicates, args.seed
    )
    cliff = build_cliff_table(ogg, args.bootstrap_replicates, args.seed + 1)
    genome, pools = pd.DataFrame(genome_rows).sort_values("Species"), pd.DataFrame(pool_rows).sort_values("Species")
    outputs = {
        "matched_background_overall.tsv": overall,
        "matched_background_by_species.tsv": by_species,
        "matched_background_replicates.tsv": replicates,
        "matched_background_representative_mapping.tsv": first_mapping,
        "background_control_pool_by_species.tsv": pools,
        "per_genome_TE_annotation_summary.tsv": genome,
        "cliffs_delta_bootstrap_ci.tsv": cliff,
    }
    for name, frame in outputs.items():
        frame.to_csv(args.output_root / name, sep="\t", index=False)
    summary = {
        "bZIP_genes": len(all_bzip),
        "bZIP_TE_positive": int(all_bzip["TE_positive"].sum()),
        "bZIP_TE_positive_fraction": float(all_bzip["TE_positive"].mean()),
        "eligible_non_bZIP_controls": len(all_controls),
        "matching_replicates": args.matching_replicates,
        "bootstrap_replicates": args.bootstrap_replicates,
        "seed": args.seed,
        "species": len(genome),
        "per_replicate_control_n": int(replicates["matched_control_genes_n"].iloc[0]),
        "matched_control_fraction_mean": float(replicates["TE_positive_control_fraction"].mean()),
        "matched_control_fraction_ci": [
            float(replicates["TE_positive_control_fraction"].quantile(0.025)),
            float(replicates["TE_positive_control_fraction"].quantile(0.975)),
        ],
        "all_passed": bool(
            len(genome) == 19
            and len(all_bzip) == 1764
            and int(all_bzip["TE_positive"].sum()) == 1274
            and replicates["matched_control_genes_n"].nunique() == 1
            and int(replicates["matched_control_genes_n"].iloc[0]) == 1764
            and len(cliff) == 3
            and (genome["nonredundant_TE_covered_bp"] <= genome["assembly_bp"]).all()
        ),
    }
    (args.output_root / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if not summary["all_passed"]:
        raise ValueError(f"QA failed: {summary}")
    print(json.dumps(summary, indent=2))

if __name__ == "__main__":
    main()
