#!/usr/bin/env python3
"""Direct alignment-based audit of group D bZIP heptad periodicity.

This analysis is deliberately independent of treating a MEME motif call as
evidence for, or against, a leucine zipper.  It uses the PF00170.27 profile to
align the bZIP domain, scores the canonical heptad a/d positions in HMM match
coordinates, and then stratifies the direct alignment result by revised MEME
Motif 2 status.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import fisher_exact, wilcoxon


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["font.size"] = 7
plt.rcParams["axes.linewidth"] = 0.8
plt.rcParams["axes.spines.top"] = False
plt.rcParams["axes.spines.right"] = False
plt.rcParams["legend.frameon"] = False


HMM_LENGTH = 64
# PF00170.27 consensus: d-position leucines at 31, 38, 45, 52, and 59.
D_STATES = (31, 38, 45, 52, 59)
# The complementary hydrophobic a phase in the same coiled-coil region.
A_STATES = (28, 35, 42, 49, 56, 63)
ZIPPER_STATES = tuple(range(28, 64))
NONCORE_STATES = tuple(
    state for state in ZIPPER_STATES if state not in set(D_STATES + A_STATES)
)
HYDROPHOBIC = frozenset("AILMFWVY")
D_OFFSETS = (0, 7, 14, 21, 28)
A_OFFSETS = (-3, 4, 11, 18, 25, 32)
WINDOW_MIN_OFFSET = -6
WINDOW_MAX_OFFSET = 34
NONCORE_OFFSETS = tuple(
    offset
    for offset in range(WINDOW_MIN_OFFSET, WINDOW_MAX_OFFSET + 1)
    if offset not in set(D_OFFSETS + A_OFFSETS)
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--primary-fasta", type=Path, required=True)
    parser.add_argument("--mast-hits", type=Path, required=True)
    parser.add_argument("--representative-mapping", type=Path, required=True)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--hmm", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def read_fasta(path: Path) -> OrderedDict[str, str]:
    records: OrderedDict[str, str] = OrderedDict()
    name: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks).upper()
                name = line[1:].split()[0]
                if name in records:
                    raise ValueError(f"Duplicate FASTA identifier: {name}")
                chunks = []
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def write_fasta(path: Path, records: OrderedDict[str, str], width: int = 80) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def run_command(command: list[str], stdout_path: Path | None = None) -> None:
    if stdout_path is None:
        completed = subprocess.run(command, check=False, text=True, capture_output=True)
        if completed.returncode != 0:
            raise RuntimeError(
                f"Command failed ({completed.returncode}): {' '.join(command)}\n"
                f"STDOUT:\n{completed.stdout}\nSTDERR:\n{completed.stderr}"
            )
    else:
        with stdout_path.open("w", encoding="utf-8", newline="\n") as output:
            completed = subprocess.run(command, check=False, text=True, stdout=output, stderr=subprocess.PIPE)
        if completed.returncode != 0:
            raise RuntimeError(
                f"Command failed ({completed.returncode}): {' '.join(command)}\n"
                f"STDERR:\n{completed.stderr}"
            )


def parse_stockholm(path: Path) -> tuple[OrderedDict[str, str], str]:
    sequences: OrderedDict[str, list[str]] = OrderedDict()
    rf_chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line == "//":
                continue
            if line.startswith("#=GC RF"):
                rf_chunks.append(line.split(maxsplit=2)[2])
                continue
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                sequences.setdefault(parts[0], []).append(parts[1])
    joined = OrderedDict((name, "".join(chunks)) for name, chunks in sequences.items())
    rf = "".join(rf_chunks)
    if not joined or not rf:
        raise ValueError(f"Could not parse Stockholm alignment: {path}")
    lengths = {len(sequence) for sequence in joined.values()}
    if lengths != {len(rf)}:
        raise ValueError(f"Stockholm alignment/RF lengths disagree: {lengths}, RF={len(rf)}")
    return joined, rf


def parse_domtblout(path: Path) -> dict[str, dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.split(maxsplit=22)
            if len(fields) < 22:
                continue
            row = {
                "sequence_id": fields[0],
                "sequence_length": int(fields[2]),
                "full_evalue": float(fields[6]),
                "full_score": float(fields[7]),
                "domain_i_evalue": float(fields[12]),
                "domain_score": float(fields[13]),
                "hmm_from": int(fields[15]),
                "hmm_to": int(fields[16]),
                "ali_from": int(fields[17]),
                "ali_to": int(fields[18]),
                "env_from": int(fields[19]),
                "env_to": int(fields[20]),
                "accuracy": float(fields[21]),
            }
            prior = best.get(fields[0])
            if prior is None or float(row["domain_score"]) > float(prior["domain_score"]):
                best[fields[0]] = row
    return best


def alignment_maps(
    aligned: OrderedDict[str, str], rf: str
) -> tuple[dict[str, dict[int, str]], dict[str, dict[int, int]], OrderedDict[str, str]]:
    match_columns = [index for index, marker in enumerate(rf) if marker not in ".-~"]
    if len(match_columns) != HMM_LENGTH:
        raise ValueError(f"Expected {HMM_LENGTH} PF00170 match states, found {len(match_columns)}")
    state_residues: dict[str, dict[int, str]] = {}
    residue_to_state: dict[str, dict[int, int]] = {}
    match_fasta: OrderedDict[str, str] = OrderedDict()
    for name, sequence in aligned.items():
        residues: dict[int, str] = {}
        residue_map: dict[int, int] = {}
        residue_number = 0
        state_number_by_column = {column: state for state, column in enumerate(match_columns, start=1)}
        for column, char in enumerate(sequence):
            if char.isalpha():
                residue_number += 1
            if column in state_number_by_column:
                state = state_number_by_column[column]
                aa = char.upper() if char.isalpha() else "-"
                residues[state] = aa
                if char.isalpha():
                    residue_map[residue_number] = state
        state_residues[name] = residues
        residue_to_state[name] = residue_map
        match_fasta[name] = "".join(residues[state] for state in range(1, HMM_LENGTH + 1))
    return state_residues, residue_to_state, match_fasta


def classify_support(d_available: int, d_hydrophobic: int) -> str:
    if d_available >= 4 and d_hydrophobic >= 4:
        return "strong"
    if d_available >= 3 and d_hydrophobic >= 3:
        return "moderate"
    return "limited"


def hydrophobic_fraction(residues: dict[int, str], states: Iterable[int]) -> tuple[int, int, float]:
    observed = [residues[state] for state in states if residues.get(state, "-") != "-"]
    hydrophobic = sum(aa in HYDROPHOBIC for aa in observed)
    fraction = hydrophobic / len(observed) if observed else math.nan
    return len(observed), hydrophobic, fraction


def format_float(value: float, digits: int = 6) -> str:
    return "" if math.isnan(value) else f"{value:.{digits}f}"


def bootstrap_mean_ci(values: list[float], seed: int = 20260808) -> tuple[float, float, float]:
    array = np.asarray([value for value in values if not math.isnan(value)], dtype=float)
    if array.size == 0:
        return math.nan, math.nan, math.nan
    rng = np.random.default_rng(seed)
    boot = rng.choice(array, size=(10000, array.size), replace=True).mean(axis=1)
    return float(array.mean()), float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.03, label, transform=ax.transAxes, fontsize=8, fontweight="bold", va="bottom")


def ptbzip_number(value: str) -> int:
    match = re.search(r"(\d+)$", value)
    if match is None:
        raise ValueError(f"PtbZIP identifier has no numeric suffix: {value}")
    return int(match.group(1))


def save_figure(fig: plt.Figure, base: Path) -> None:
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".png"), dpi=400, bbox_inches="tight")
    fig.savefig(base.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)


def make_summary_figure(
    output_dir: Path,
    per_protein: list[dict[str, object]],
    pangene_rows: list[dict[str, object]],
    motif2_span_summary: dict[str, object],
) -> None:
    fig = plt.figure(figsize=(7.2, 5.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[0.8, 3.2], width_ratios=[1.45, 1.0])

    ax_a = fig.add_subplot(grid[0, :])
    ax_a.set_xlim(0, 65)
    ax_a.set_ylim(0, 1)
    ax_a.axis("off")
    ax_a.add_patch(plt.Rectangle((1, 0.34), 27, 0.28, color="#CFE7F1", ec="none"))
    ax_a.add_patch(plt.Rectangle((28, 0.34), 36, 0.28, color="#F4C3A2", ec="none"))
    ax_a.text(14.5, 0.48, "Basic region", ha="center", va="center", fontweight="bold")
    ax_a.text(46, 0.48, "Leucine zipper", ha="center", va="center", fontweight="bold")
    for index, state in enumerate(D_STATES, start=1):
        ax_a.plot([state, state], [0.26, 0.72], color="#B64342", lw=1.3)
        ax_a.text(state, 0.76, f"d{index}", ha="center", va="bottom", color="#B64342", fontsize=6.5)
    for state in A_STATES:
        ax_a.plot([state, state], [0.34, 0.62], color="#3775BA", lw=0.6, alpha=0.8)
    mapped = motif2_span_summary.get("median_hmm_span")
    if isinstance(mapped, list) and len(mapped) == 2:
        start, end = mapped
        ax_a.annotate(
            "Motif 2 median aligned span",
            xy=((start + end) / 2, 0.19),
            xytext=((start + end) / 2, 0.03),
            ha="center",
            va="top",
            arrowprops={"arrowstyle": "-[,widthB=3.5,lengthB=0.5", "lw": 0.8, "color": "#0F4D92"},
            color="#0F4D92",
            fontsize=6.5,
        )
    ax_a.text(1, 0.9, "PF00170.27 match-state coordinate", ha="left", va="center", fontsize=6.5)
    add_panel_label(ax_a, "a")

    ax_b = fig.add_subplot(grid[1, 0])
    labels = [str(row["group_label"]) for row in pangene_rows]
    matrix = np.array(
        [[float(row[f"d{index}_hydrophobic_percent"]) for index in range(1, 6)] + [float(row["motif2_percent"])] for row in pangene_rows]
    )
    cmap = LinearSegmentedColormap.from_list("heptad", ["#F5F5F5", "#B4CBE6", "#0F4D92"])
    cmap.set_bad("white")
    image = ax_b.imshow(np.ma.masked_invalid(matrix), vmin=0, vmax=100, cmap=cmap, aspect="auto")
    ax_b.set_xticks(range(6), ["d1", "d2", "d3", "d4", "d5", "Motif 2"])
    ax_b.set_yticks(range(len(labels)), labels)
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            if math.isnan(value):
                ax_b.text(column_index, row_index, "NA", ha="center", va="center", fontsize=5.2, color="#767676")
                continue
            color = "white" if value >= 62 else "#272727"
            ax_b.text(column_index, row_index, f"{value:.0f}", ha="center", va="center", fontsize=5.2, color=color)
    ax_b.axvline(4.5, color="white", lw=2)
    colorbar = fig.colorbar(image, ax=ax_b, fraction=0.035, pad=0.02)
    colorbar.set_label("Hydrophobic residues or Motif 2-positive proteins (%)")
    ax_b.set_xlabel("Fixed heptad d positions in the PF00170 alignment")
    add_panel_label(ax_b, "b")

    ax_c = fig.add_subplot(grid[1, 1])
    groups = [
        ("Populus\nMotif 2+", "Populus_Motif2_positive", "#3775BA"),
        ("Populus\nMotif 2−", "Populus_Motif2_negative", "#8F8F8F"),
        ("Arabidopsis\ngroup D", "Arabidopsis_D", "#D98A82"),
        ("Rice\ngroup D", "Rice_D", "#85AA8B"),
    ]
    rng = np.random.default_rng(20260808)
    box_data: list[list[float]] = []
    for _, group_key, _ in groups:
        values = [
            float(row["periodicity_score"])
            for row in per_protein
            if row["analysis_group"] == group_key and row["periodicity_score"] != ""
        ]
        box_data.append(values)
    box = ax_c.boxplot(box_data, widths=0.55, patch_artist=True, showfliers=False, medianprops={"color": "#272727", "lw": 1})
    for patch, (_, _, color) in zip(box["boxes"], groups):
        patch.set_facecolor(color)
        patch.set_alpha(0.55)
        patch.set_edgecolor("#4D4D4D")
    for index, (values, (_, _, color)) in enumerate(zip(box_data, groups), start=1):
        x = rng.normal(index, 0.055, len(values))
        ax_c.scatter(x, values, s=6, color=color, alpha=0.45, edgecolors="none")
        ax_c.text(index, -0.47, f"n={len(values)}", ha="center", va="top", fontsize=5.8)
    ax_c.axhline(0, color="#767676", lw=0.8, ls="--")
    ax_c.set_xticks(range(1, 5), [group[0] for group in groups])
    ax_c.set_ylabel("Heptad periodicity score\n(d-position hydrophobic fraction − non-core fraction)")
    ax_c.set_ylim(-0.45, 1.05)
    add_panel_label(ax_c, "c")

    save_figure(fig, output_dir / "Figure_group_D_heptad_periodicity_summary")


def make_alignment_figure(
    output_dir: Path,
    per_protein: list[dict[str, object]],
) -> None:
    representatives = [row for row in per_protein if row["origin"] == "Populus" and row["is_representative"] == "yes"]
    representatives.sort(key=lambda row: ptbzip_number(str(row["PtbZIP_ID"])))
    arabidopsis = [row for row in per_protein if row["origin"] == "Arabidopsis"]
    arabidopsis.sort(key=lambda row: str(row["sequence_id"]))
    display_rows = representatives + arabidopsis
    offsets = list(range(WINDOW_MIN_OFFSET, WINDOW_MAX_OFFSET + 1))
    fig_height = max(5.5, 0.255 * len(display_rows) + 1.6)
    fig, ax = plt.subplots(figsize=(7.2, fig_height))
    ax.set_xlim(-0.5, len(offsets) - 0.5)
    ax.set_ylim(len(display_rows) - 0.5, -0.5)
    tick_columns = [index for index, offset in enumerate(offsets) if offset % 7 == 0 or offset in {WINDOW_MIN_OFFSET, WINDOW_MAX_OFFSET}]
    ax.set_xticks(tick_columns, [offsets[index] for index in tick_columns], fontsize=5.3)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=2, pad=1)
    labels = [str(row["PtbZIP_ID"] if row["origin"] == "Populus" else row["sequence_id"]) for row in display_rows]
    ax.set_yticks(range(len(display_rows)), labels, fontsize=6)
    for column, offset in enumerate(offsets):
        if offset in D_OFFSETS:
            ax.axvspan(column - 0.5, column + 0.5, color="#F4C3A2", alpha=0.65, zorder=0)
        elif offset in A_OFFSETS:
            ax.axvspan(column - 0.5, column + 0.5, color="#CFE7F1", alpha=0.6, zorder=0)
    for row_index, row in enumerate(display_rows):
        window = str(row["zipper_window_sequence"])
        for column, (offset, aa) in enumerate(zip(offsets, window)):
            color = "#B64342" if offset in D_OFFSETS and aa in HYDROPHOBIC else "#272727"
            weight = "bold" if offset in D_OFFSETS else "normal"
            ax.text(column, row_index, aa, ha="center", va="center", fontsize=5.4, color=color, fontweight=weight)
        if row["origin"] == "Populus" and row["motif2_present"] == "yes":
            ax.scatter(len(offsets) - 0.05, row_index, s=12, color="#3775BA", clip_on=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Residue offset from PF00170-anchored d1 (a positions blue; d positions orange)")
    ax.xaxis.set_label_position("top")
    ax.text(1.01, 1.025, "● Motif 2 detected", transform=ax.transAxes, color="#3775BA", fontsize=6, ha="left", va="bottom")
    ax.text(-0.19, 1.025, "Group D representatives", transform=ax.transAxes, fontsize=8, fontweight="bold", ha="left", va="bottom")
    save_figure(fig, output_dir / "Figure_group_D_representative_heptad_alignment")


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(f"Output directory already exists; choose a new root: {args.output_dir}")
    args.output_dir.mkdir(parents=True)
    source_snapshot = args.output_dir / "source_snapshot"
    source_snapshot.mkdir()

    mapping_rows = read_tsv(args.mapping)
    primary_sequences = read_fasta(args.primary_fasta)
    reference_sequences = read_fasta(args.reference_fasta)
    mast_hits = read_tsv(args.mast_hits)
    representative_rows = read_tsv(args.representative_mapping)

    group_d_mapping = [
        row
        for row in mapping_rows
        if row["subfamily"] == "D" and row["assignment_status"] == "assigned" and row["included_primary_1764"] == "yes"
    ]
    if len(group_d_mapping) != 255:
        raise ValueError(f"Expected 255 assigned primary group D proteins, found {len(group_d_mapping)}")
    missing = [row["upload_id"] for row in group_d_mapping if row["upload_id"] not in primary_sequences]
    if missing:
        raise ValueError(f"Group D identifiers missing from primary FASTA: {missing[:10]}")

    group_d_references = OrderedDict((name, sequence) for name, sequence in reference_sequences.items() if name.endswith("_D"))
    at_count = sum(name.startswith("At") for name in group_d_references)
    os_count = sum(name.startswith("Os") for name in group_d_references)
    if (at_count, os_count) != (10, 16):
        raise ValueError(f"Expected 10 Arabidopsis and 16 rice group D references, found {at_count} and {os_count}")

    motif2_hits = {
        row["upload_id"]: row
        for row in mast_hits
        if row["subfamily"] == "D" and row["motif_number"] == "2"
    }
    if len(motif2_hits) != 37:
        raise ValueError(f"Expected 37 group D Motif 2 hits, found {len(motif2_hits)}")

    representative_by_member = {
        row["representative_member_id"]: row
        for row in representative_rows
        if row["subfamily"] == "D"
    }
    if len(representative_by_member) != 12:
        raise ValueError(f"Expected 12 official group D representatives, found {len(representative_by_member)}")

    combined: OrderedDict[str, str] = OrderedDict()
    metadata: dict[str, dict[str, object]] = {}
    for row in sorted(group_d_mapping, key=lambda item: item["upload_id"]):
        sequence_id = row["upload_id"]
        combined[sequence_id] = primary_sequences[sequence_id]
        hit = motif2_hits.get(sequence_id)
        metadata[sequence_id] = {
            "sequence_id": sequence_id,
            "origin": "Populus",
            "species": row["species"],
            "original_id": row["original_id"],
            "PtbZIP_ID": row["PtbZIP_ID"],
            "pan_genome_class": row["pan_genome_class"],
            "sequence_length": len(primary_sequences[sequence_id]),
            "motif2_present": "yes" if hit else "no",
            "motif2_hit_position_1based": int(hit["hit_position_1based"]) if hit else "",
            "motif2_hit_pvalue": hit["hit_pvalue"] if hit else "",
            "is_representative": "yes" if row["original_id"] in representative_by_member else "no",
        }
    for sequence_id, sequence in group_d_references.items():
        combined[sequence_id] = sequence
        origin = "Arabidopsis" if sequence_id.startswith("At") else "Rice"
        metadata[sequence_id] = {
            "sequence_id": sequence_id,
            "origin": origin,
            "species": "Arabidopsis_thaliana" if origin == "Arabidopsis" else "Oryza_sativa",
            "original_id": sequence_id,
            "PtbZIP_ID": "",
            "pan_genome_class": "reference",
            "sequence_length": len(sequence),
            "motif2_present": "not_scanned",
            "motif2_hit_position_1based": "",
            "motif2_hit_pvalue": "",
            "is_representative": "reference",
        }

    combined_fasta = args.output_dir / "group_D_Populus255_plus_AtOs26.full_length.fasta"
    write_fasta(combined_fasta, combined)
    shutil.copy2(args.hmm, source_snapshot / "PF00170.27.hmm")
    shutil.copy2(args.reference_fasta, source_snapshot / "AtOs_bZIP_group_reference_source.fasta")

    domtblout = args.output_dir / "PF00170_group_D.domtblout"
    hmmsearch_text = args.output_dir / "PF00170_group_D.hmmsearch.txt"
    run_command(
        [
            "hmmsearch",
            "--cpu",
            str(args.threads),
            "--domtblout",
            str(domtblout),
            "-o",
            str(hmmsearch_text),
            str(args.hmm),
            str(combined_fasta),
        ]
    )
    stockholm = args.output_dir / "PF00170_group_D.hmmalign.sto"
    run_command(
        ["hmmalign", "--outformat", "Stockholm", str(args.hmm), str(combined_fasta)],
        stdout_path=stockholm,
    )

    aligned, rf = parse_stockholm(stockholm)
    if set(aligned) != set(combined):
        raise ValueError(f"Alignment identifiers differ from combined FASTA: aligned={len(aligned)}, input={len(combined)}")
    state_residues, residue_to_state, match_fasta = alignment_maps(aligned, rf)
    match_fasta_path = args.output_dir / "PF00170_group_D.match_states_64.fasta"
    write_fasta(match_fasta_path, match_fasta)
    domains = parse_domtblout(domtblout)

    motif_domain_rows: list[dict[str, object]] = []
    for hit in mast_hits:
        if hit["subfamily"] != "D":
            continue
        sequence_id = hit["upload_id"]
        domain = domains.get(sequence_id, {})
        start = int(hit["hit_position_1based"])
        end = start + int(hit["motif_width"]) - 1
        ali_from = int(domain.get("ali_from", 0))
        ali_to = int(domain.get("ali_to", 0))
        env_from = int(domain.get("env_from", 0))
        env_to = int(domain.get("env_to", 0))
        overlap_ali = max(0, min(end, ali_to) - max(start, ali_from) + 1) if ali_from and ali_to else 0
        overlap_env = max(0, min(end, env_to) - max(start, env_from) + 1) if env_from and env_to else 0
        motif_domain_rows.append(
            {
                "sequence_id": sequence_id,
                "PtbZIP_ID": metadata[sequence_id]["PtbZIP_ID"],
                "motif_number": int(hit["motif_number"]),
                "motif_width": int(hit["motif_width"]),
                "hit_start_1based": start,
                "hit_end_1based": end,
                "PF00170_ali_from": ali_from or "",
                "PF00170_ali_to": ali_to or "",
                "PF00170_env_from": env_from or "",
                "PF00170_env_to": env_to or "",
                "overlap_with_alignment_residues": overlap_ali,
                "overlap_with_envelope_residues": overlap_env,
            }
        )
    write_tsv(
        args.output_dir / "group_D_all_motif_hits_PF00170_overlap.tsv",
        motif_domain_rows,
        list(motif_domain_rows[0]),
    )
    motif_overlap_summary: list[dict[str, object]] = []
    for motif_number in range(1, 11):
        selected = [row for row in motif_domain_rows if row["motif_number"] == motif_number]
        motif_overlap_summary.append(
            {
                "motif_number": motif_number,
                "group_D_genes_with_hit": len({row["sequence_id"] for row in selected}),
                "hits_total": len(selected),
                "hits_overlapping_PF00170_alignment": sum(int(row["overlap_with_alignment_residues"]) > 0 for row in selected),
                "hits_overlapping_PF00170_envelope": sum(int(row["overlap_with_envelope_residues"]) > 0 for row in selected),
                "genes_with_alignment_overlapping_hit": len({row["sequence_id"] for row in selected if int(row["overlap_with_alignment_residues"]) > 0}),
                "genes_with_envelope_overlapping_hit": len({row["sequence_id"] for row in selected if int(row["overlap_with_envelope_residues"]) > 0}),
            }
        )
    write_tsv(
        args.output_dir / "group_D_motif_PF00170_overlap_summary.tsv",
        motif_overlap_summary,
        list(motif_overlap_summary[0]),
    )

    motif2_mapping_rows: list[dict[str, object]] = []
    motif2_spans: list[tuple[int, int]] = []
    motif2_l_state_triplets: list[tuple[int | None, int | None, int | None]] = []
    for sequence_id, hit in motif2_hits.items():
        start = int(hit["hit_position_1based"])
        end = start + int(hit["motif_width"]) - 1
        domain = domains.get(sequence_id, {})
        ali_from = int(domain.get("ali_from", 0))
        ali_to = int(domain.get("ali_to", 0))
        env_from = int(domain.get("env_from", 0))
        env_to = int(domain.get("env_to", 0))
        overlap_ali = max(0, min(end, ali_to) - max(start, ali_from) + 1) if ali_from and ali_to else 0
        overlap_env = max(0, min(end, env_to) - max(start, env_from) + 1) if env_from and env_to else 0
        mapped_states = [
            residue_to_state[sequence_id][position]
            for position in range(start, end + 1)
            if position in residue_to_state[sequence_id]
        ]
        if mapped_states:
            span = (min(mapped_states), max(mapped_states))
            motif2_spans.append(span)
        else:
            span = (None, None)
        l_positions = (start + 4, start + 11, start + 18)
        l_states = tuple(residue_to_state[sequence_id].get(position) for position in l_positions)
        motif2_l_state_triplets.append(l_states)
        motif2_mapping_rows.append(
            {
                "sequence_id": sequence_id,
                "PtbZIP_ID": metadata[sequence_id]["PtbZIP_ID"],
                "motif2_start_1based": start,
                "motif2_end_1based": end,
                "PF00170_ali_from": ali_from or "",
                "PF00170_ali_to": ali_to or "",
                "PF00170_env_from": env_from or "",
                "PF00170_env_to": env_to or "",
                "overlap_with_PF00170_alignment_residues": overlap_ali,
                "overlaps_PF00170_alignment": "yes" if overlap_ali else "no",
                "overlap_with_PF00170_envelope_residues": overlap_env,
                "overlaps_PF00170_envelope": "yes" if overlap_env else "no",
                "mapped_PF00170_match_states_count": len(mapped_states),
                "hmm_span_start": span[0] or "",
                "hmm_span_end": span[1] or "",
                "motif2_consensus_L1_hmm_state": l_states[0] or "",
                "motif2_consensus_L2_hmm_state": l_states[1] or "",
                "motif2_consensus_L3_hmm_state": l_states[2] or "",
            }
        )
    write_tsv(
        args.output_dir / "motif2_to_PF00170_alignment_mapping.tsv",
        motif2_mapping_rows,
        list(motif2_mapping_rows[0]),
    )
    median_span = [int(np.median([span[index] for span in motif2_spans])) for index in (0, 1)] if motif2_spans else None
    l_state_modes: list[int | None] = []
    for index in range(3):
        counts: dict[int, int] = defaultdict(int)
        for triplet in motif2_l_state_triplets:
            if triplet[index] is not None:
                counts[int(triplet[index])] += 1
        l_state_modes.append(max(counts, key=counts.get) if counts else None)
    motif2_span_summary = {
        "motif2_positive_group_D_proteins": len(motif2_hits),
        "motif2_hits_overlapping_PF00170_alignment": sum(row["overlaps_PF00170_alignment"] == "yes" for row in motif2_mapping_rows),
        "motif2_hits_overlapping_PF00170_envelope": sum(row["overlaps_PF00170_envelope"] == "yes" for row in motif2_mapping_rows),
        "median_hmm_span": median_span,
        "modal_hmm_states_for_motif2_consensus_L_positions_5_12_19": l_state_modes,
        "canonical_d_states": list(D_STATES),
    }
    (args.output_dir / "motif2_to_PF00170_alignment_summary.json").write_text(
        json.dumps(motif2_span_summary, indent=2) + "\n", encoding="utf-8"
    )

    per_protein: list[dict[str, object]] = []
    for sequence_id in combined:
        meta = metadata[sequence_id]
        residues = state_residues[sequence_id]
        anchor_positions = [position for position, state in residue_to_state[sequence_id].items() if state == D_STATES[0]]
        anchor_position = anchor_positions[0] if len(anchor_positions) == 1 else None
        full_sequence = combined[sequence_id]
        extended: dict[int, str] = {}
        for offset in range(WINDOW_MIN_OFFSET, WINDOW_MAX_OFFSET + 1):
            raw_index = (anchor_position - 1 + offset) if anchor_position is not None else -1
            extended[offset] = full_sequence[raw_index] if 0 <= raw_index < len(full_sequence) else "-"
        d_available, d_hydrophobic, d_fraction = hydrophobic_fraction(extended, D_OFFSETS)
        a_available, a_hydrophobic, a_fraction = hydrophobic_fraction(extended, A_OFFSETS)
        noncore_available, noncore_hydrophobic, noncore_fraction = hydrophobic_fraction(extended, NONCORE_OFFSETS)
        periodicity_score = d_fraction - noncore_fraction if not math.isnan(d_fraction) and not math.isnan(noncore_fraction) else math.nan
        if meta["origin"] == "Populus":
            analysis_group = "Populus_Motif2_positive" if meta["motif2_present"] == "yes" else "Populus_Motif2_negative"
        else:
            analysis_group = f"{meta['origin']}_D"
        domain = domains.get(sequence_id, {})
        row: dict[str, object] = {
            **meta,
            "analysis_group": analysis_group,
            "d1_anchor_full_sequence_position_1based": anchor_position or "",
            "zipper_window_start_offset": WINDOW_MIN_OFFSET,
            "zipper_window_end_offset": WINDOW_MAX_OFFSET,
            "zipper_window_sequence": "".join(extended[offset] for offset in range(WINDOW_MIN_OFFSET, WINDOW_MAX_OFFSET + 1)),
            "domain_i_evalue": domain.get("domain_i_evalue", ""),
            "domain_score": domain.get("domain_score", ""),
            "hmm_from": domain.get("hmm_from", ""),
            "hmm_to": domain.get("hmm_to", ""),
            "ali_from": domain.get("ali_from", ""),
            "ali_to": domain.get("ali_to", ""),
            "d_available": d_available,
            "d_hydrophobic_count": d_hydrophobic,
            "d_leucine_count": sum(residues[state] == "L" for state in D_STATES),
            "d_hydrophobic_fraction": format_float(d_fraction),
            "a_available": a_available,
            "a_hydrophobic_count": a_hydrophobic,
            "a_hydrophobic_fraction": format_float(a_fraction),
            "noncore_available": noncore_available,
            "noncore_hydrophobic_count": noncore_hydrophobic,
            "noncore_hydrophobic_fraction": format_float(noncore_fraction),
            "periodicity_score": format_float(periodicity_score),
            "heptad_support": classify_support(d_available, d_hydrophobic),
        }
        for index, (state, offset) in enumerate(zip(D_STATES, D_OFFSETS), start=1):
            aa = extended[offset]
            row[f"d{index}_state"] = state
            row[f"d{index}_offset_from_anchor"] = offset
            row[f"d{index}_residue"] = aa
            row[f"d{index}_hydrophobic"] = "yes" if aa in HYDROPHOBIC else "no"
        per_protein.append(row)

    per_protein_fields = list(per_protein[0])
    write_tsv(args.output_dir / "group_D_heptad_per_protein.tsv", per_protein, per_protein_fields)

    summary_groups = [
        ("Populus_all_D", lambda row: row["origin"] == "Populus"),
        ("Populus_Motif2_positive", lambda row: row["analysis_group"] == "Populus_Motif2_positive"),
        ("Populus_Motif2_negative", lambda row: row["analysis_group"] == "Populus_Motif2_negative"),
        ("Arabidopsis_D", lambda row: row["analysis_group"] == "Arabidopsis_D"),
        ("Rice_D", lambda row: row["analysis_group"] == "Rice_D"),
        ("AtOs_D_all", lambda row: row["origin"] in {"Arabidopsis", "Rice"}),
    ]
    summary_rows: list[dict[str, object]] = []
    for group_name, predicate in summary_groups:
        selected = [row for row in per_protein if predicate(row)]
        d_values = [float(row["d_hydrophobic_fraction"]) for row in selected if row["d_hydrophobic_fraction"] != ""]
        score_values = [float(row["periodicity_score"]) for row in selected if row["periodicity_score"] != ""]
        score_mean, score_low, score_high = bootstrap_mean_ci(score_values)
        paired = [
            (float(row["d_hydrophobic_fraction"]), float(row["noncore_hydrophobic_fraction"]))
            for row in selected
            if row["d_hydrophobic_fraction"] != "" and row["noncore_hydrophobic_fraction"] != ""
        ]
        if paired and any(abs(d - n) > 1e-12 for d, n in paired):
            wilcoxon_result = wilcoxon([d for d, _ in paired], [n for _, n in paired], alternative="greater")
            wilcoxon_stat = float(wilcoxon_result.statistic)
            wilcoxon_p = float(wilcoxon_result.pvalue)
        else:
            wilcoxon_stat = math.nan
            wilcoxon_p = math.nan
        summary_rows.append(
            {
                "analysis_group": group_name,
                "n": len(selected),
                "motif2_positive_n": sum(row["motif2_present"] == "yes" for row in selected),
                "motif2_positive_percent": format_float(100 * sum(row["motif2_present"] == "yes" for row in selected) / len(selected)) if selected and all(row["origin"] == "Populus" for row in selected) else "",
                "heptad_evaluable_n": sum(row["d1_anchor_full_sequence_position_1based"] != "" for row in selected),
                "mean_d_hydrophobic_fraction": format_float(float(np.mean(d_values)) if d_values else math.nan),
                "proteins_with_at_least_3_of_5_hydrophobic_d_positions": sum(int(row["d_hydrophobic_count"]) >= 3 for row in selected),
                "alignment_supported_at_least3_percent_of_evaluable": format_float(100 * sum(int(row["d_hydrophobic_count"]) >= 3 for row in selected) / sum(row["d1_anchor_full_sequence_position_1based"] != "" for row in selected)) if any(row["d1_anchor_full_sequence_position_1based"] != "" for row in selected) else "",
                "proteins_with_at_least_4_of_5_hydrophobic_d_positions": sum(int(row["d_hydrophobic_count"]) >= 4 for row in selected),
                "proteins_with_5_of_5_hydrophobic_d_positions": sum(int(row["d_hydrophobic_count"]) == 5 for row in selected),
                "strong_support_n": sum(row["heptad_support"] == "strong" for row in selected),
                "moderate_support_n": sum(row["heptad_support"] == "moderate" for row in selected),
                "limited_support_n": sum(row["heptad_support"] == "limited" for row in selected),
                "mean_periodicity_score": format_float(score_mean),
                "periodicity_score_bootstrap_95CI_low": format_float(score_low),
                "periodicity_score_bootstrap_95CI_high": format_float(score_high),
                "wilcoxon_d_vs_noncore_statistic": format_float(wilcoxon_stat),
                "wilcoxon_d_vs_noncore_one_sided_pvalue": f"{wilcoxon_p:.6g}" if not math.isnan(wilcoxon_p) else "",
            }
        )
    write_tsv(args.output_dir / "group_D_heptad_group_summary.tsv", summary_rows, list(summary_rows[0]))

    motif2_positive = [row for row in per_protein if row["analysis_group"] == "Populus_Motif2_positive" and row["d1_anchor_full_sequence_position_1based"] != ""]
    motif2_negative = [row for row in per_protein if row["analysis_group"] == "Populus_Motif2_negative" and row["d1_anchor_full_sequence_position_1based"] != ""]
    contingency = [
        [sum(int(row["d_hydrophobic_count"]) >= 3 for row in motif2_positive), sum(int(row["d_hydrophobic_count"]) < 3 for row in motif2_positive)],
        [sum(int(row["d_hydrophobic_count"]) >= 3 for row in motif2_negative), sum(int(row["d_hydrophobic_count"]) < 3 for row in motif2_negative)],
    ]
    odds_ratio, fisher_p = fisher_exact(contingency, alternative="two-sided")
    contingency_row = {
        "motif2_positive_evaluable": len(motif2_positive),
        "motif2_positive_alignment_supported_at_least3": contingency[0][0],
        "motif2_positive_not_supported": contingency[0][1],
        "motif2_negative_evaluable": len(motif2_negative),
        "motif2_negative_alignment_supported_at_least3": contingency[1][0],
        "motif2_negative_not_supported": contingency[1][1],
        "fisher_exact_odds_ratio": odds_ratio,
        "fisher_exact_two_sided_pvalue": fisher_p,
    }
    write_tsv(args.output_dir / "motif2_vs_direct_heptad_support.tsv", [contingency_row], list(contingency_row))

    pangene_rows: list[dict[str, object]] = []
    labels_and_predicates: list[tuple[str, object]] = []
    for pangene in sorted({str(row["PtbZIP_ID"]) for row in per_protein if row["origin"] == "Populus"}, key=ptbzip_number):
        labels_and_predicates.append((pangene, lambda row, pangene=pangene: row["PtbZIP_ID"] == pangene))
    labels_and_predicates.extend(
        [
            ("Arabidopsis D", lambda row: row["origin"] == "Arabidopsis"),
            ("Rice D", lambda row: row["origin"] == "Rice"),
        ]
    )
    for label, predicate in labels_and_predicates:
        selected = [row for row in per_protein if predicate(row)]
        row_out: dict[str, object] = {
            "group_label": label,
            "n": len(selected),
            "motif2_positive_n": sum(row["motif2_present"] == "yes" for row in selected),
            "motif2_percent": 100 * sum(row["motif2_present"] == "yes" for row in selected) / len(selected) if selected and label.startswith("PtbZIP") else math.nan,
        }
        for index in range(1, 6):
            observed = [row for row in selected if row[f"d{index}_residue"] != "-"]
            row_out[f"d{index}_hydrophobic_percent"] = 100 * sum(row[f"d{index}_hydrophobic"] == "yes" for row in observed) / len(observed) if observed else math.nan
        pangene_rows.append(row_out)
    write_tsv(args.output_dir / "group_D_heptad_summary_by_pangene.tsv", pangene_rows, list(pangene_rows[0]))

    make_summary_figure(args.output_dir, per_protein, pangene_rows, motif2_span_summary)
    make_alignment_figure(args.output_dir, per_protein)

    source_code_dir = args.output_dir / "source_code"
    source_code_dir.mkdir()
    shutil.copy2(Path(__file__), source_code_dir / Path(__file__).name)
    software_versions = {
        "python": sys.version.split()[0],
        "matplotlib": matplotlib.__version__,
        "numpy": np.__version__,
        "scipy": __import__("scipy").__version__,
        "hmmsearch": subprocess.run(["hmmsearch", "-h"], check=False, text=True, capture_output=True).stdout.splitlines()[1].strip(),
        "hmmalign": subprocess.run(["hmmalign", "-h"], check=False, text=True, capture_output=True).stdout.splitlines()[1].strip(),
        "profile": "PF00170.27; HMM length 64; GA 27.6",
    }
    (args.output_dir / "software_versions.json").write_text(json.dumps(software_versions, indent=2) + "\n", encoding="utf-8")

    qa = {
        "all_passed": True,
        "populus_group_D_proteins": sum(row["origin"] == "Populus" for row in per_protein),
        "arabidopsis_group_D_references": sum(row["origin"] == "Arabidopsis" for row in per_protein),
        "rice_group_D_references": sum(row["origin"] == "Rice" for row in per_protein),
        "motif2_positive_populus_group_D": len(motif2_hits),
        "alignment_sequences": len(aligned),
        "PF00170_match_states": HMM_LENGTH,
        "canonical_d_states": list(D_STATES),
        "canonical_a_states": list(A_STATES),
        "hydrophobic_residues": "".join(sorted(HYDROPHOBIC)),
        "d1_anchor_missing": sum(row["d1_anchor_full_sequence_position_1based"] == "" for row in per_protein),
        "motif2_median_hmm_span": median_span,
        "motif2_hits_overlapping_PF00170_alignment": motif2_span_summary["motif2_hits_overlapping_PF00170_alignment"],
        "motif2_hits_overlapping_PF00170_envelope": motif2_span_summary["motif2_hits_overlapping_PF00170_envelope"],
        "motif2_modal_L_states": l_state_modes,
    }
    expected = {
        "populus_group_D_proteins": 255,
        "arabidopsis_group_D_references": 10,
        "rice_group_D_references": 16,
        "motif2_positive_populus_group_D": 37,
        "alignment_sequences": 281,
        "PF00170_match_states": 64,
        "d1_anchor_missing": 5,
    }
    failures = [key for key, value in expected.items() if qa[key] != value]
    qa["failures"] = failures
    qa["all_passed"] = not failures
    (args.output_dir / "QA_REPORT.json").write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")

    methods = f"""# Direct group D heptad-periodicity assessment

## Scientific conclusion tested

Revised MEME Motif 2 has a leucine-heptad-like consensus, but its group D hits occur outside the PF00170 alignment and domain envelope. Direct PF00170 alignment is therefore used to test fixed seven-residue hydrophobic periodicity independently of the MEME call.

## Inputs and method

- Populus group D proteins: 255 OGG-assigned proteins from the revised 1,764-protein primary set.
- Positive reference set: 10 Arabidopsis and 16 rice group D bZIP proteins used in the revised phylogeny.
- Domain profile: PF00170.27 (64 match states; Pfam gathering threshold 27.6).
- Direct alignment: HMMER hmmalign with the PF00170.27 profile. Profile match states were extracted from the Stockholm RF annotation; insertions and terminal nonhomologous residues were excluded from the 64-state comparison while full-sequence residue numbering was retained for mapping Motif 2 coordinates.
- Canonical d positions: PF00170 match state 31 anchors d1; d2-d5 were read directly from each full-length protein at +7, +14, +21, and +28 residues. The corresponding PF00170 consensus states are {', '.join(map(str, D_STATES))}.
- Complementary a positions: PF00170 match states {', '.join(map(str, A_STATES))}.
- Hydrophobic residues: {''.join(sorted(HYDROPHOBIC))}.
- Descriptive strong support: at least four hydrophobic residues among at least four observed d positions; moderate support: at least three among at least three observed d positions. Full position-level results are retained so the conclusion does not depend only on these descriptive labels.
- Periodicity score: hydrophobic fraction at the five anchored d positions minus the hydrophobic fraction at non-a/non-d offsets from -6 through +34 around d1.
- Group tests: one-sided paired Wilcoxon signed-rank tests compare d-position and non-core hydrophobic fractions within sequences. Bootstrap confidence intervals use 10,000 resamples and seed 20260808.
- Motif-domain localization: every group D MAST hit for Motifs 1-10 was intersected with the best PF00170 alignment and envelope coordinates from hmmsearch. This distinguishes a motif call elsewhere in a full-length protein from a motif call within the canonical bZIP domain.

## Reviewer-risk boundary

This is a sequence-alignment assessment of heptad periodicity, not an experimental demonstration of dimerization or coiled-coil formation. Motif 2 is reported as a partial motif-level representation of the zipper, while the direct PF00170 alignment is the primary evidence for the leucine-zipper conclusion.
"""
    (args.output_dir / "README.md").write_text(methods, encoding="utf-8")

    manifest_rows = []
    for path in sorted(args.output_dir.rglob("*")):
        if path.is_file() and path.name != "SHA256SUMS.tsv":
            manifest_rows.append({"relative_path": path.relative_to(args.output_dir).as_posix(), "sha256": sha256(path), "bytes": path.stat().st_size})
    write_tsv(args.output_dir / "SHA256SUMS.tsv", manifest_rows, ["relative_path", "sha256", "bytes"])

    if not qa["all_passed"]:
        raise RuntimeError(f"QA failed: {failures}")
    print(json.dumps(qa, indent=2))


if __name__ == "__main__":
    main()
