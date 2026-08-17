#!/usr/bin/env python3
"""Recover run-level HISAT2 alignment rates from archived RNA-seq outputs."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


RUN_RE = re.compile(r"[/\\]([DES]RR\d+)(?:_[12])?\.fastq(?:\.gz)?", re.IGNORECASE)
RATE_RE = re.compile(r"([0-9]+(?:\.[0-9]+)?)% overall alignment rate")
REP_RE = re.compile(r"_rep(\d+)$", re.IGNORECASE)


def recover_run_qc(archive_root: Path) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for fastp_json in archive_root.rglob("fastp.json"):
        sample_dir = fastp_json.parent.parent
        hisat_summary = sample_dir / "qc" / "hisat2_summary.txt"
        if not hisat_summary.exists():
            continue
        payload = json.loads(fastp_json.read_text(encoding="utf-8"))
        command = str(payload.get("command", ""))
        run_match = RUN_RE.search(command)
        rate_match = RATE_RE.search(hisat_summary.read_text(encoding="utf-8"))
        rep_match = REP_RE.search(sample_dir.name)
        if not run_match or not rate_match:
            continue
        records.append(
            {
                "Run Accession": run_match.group(1).upper(),
                "archive_sample_name": sample_dir.name,
                "archive_biological_replicate": int(rep_match.group(1)) if rep_match else pd.NA,
                "Mapping rate (%)": float(rate_match.group(1)),
                "Mapping QC note": "Recovered from the archived HISAT2 summary",
                "HISAT2 summary path": str(hisat_summary),
            }
        )
    result = pd.DataFrame(records)
    if result.empty:
        raise RuntimeError(f"No recoverable HISAT2 summaries found under {archive_root}")
    if result["Run Accession"].duplicated().any():
        duplicated = result.loc[result["Run Accession"].duplicated(False), "Run Accession"].tolist()
        raise RuntimeError(f"Run accessions mapped more than once: {duplicated}")
    return result.sort_values("Run Accession").reset_index(drop=True)


def update_s13(s13_path: Path, run_qc: pd.DataFrame) -> pd.DataFrame:
    s13 = pd.read_csv(s13_path, sep="\t")
    lookup = run_qc.set_index("Run Accession")
    for row_index, run in s13["Run Accession"].astype(str).items():
        if run not in lookup.index:
            continue
        recovered = lookup.loc[run]
        s13.at[row_index, "Mapping rate (%)"] = recovered["Mapping rate (%)"]
        s13.at[row_index, "Mapping QC note"] = recovered["Mapping QC note"]
        s13.at[row_index, "Biological replicate"] = recovered["archive_biological_replicate"]
    return s13


def update_contrast_qc(qc_path: Path, s13: pd.DataFrame) -> pd.DataFrame:
    qc = pd.read_csv(qc_path, sep="\t")
    species_lookup = {
        "Populus_alba_var__pyramidalis": "Populus alba var. pyramidalis",
        "Populus_deltoides": "Populus deltoides",
        "Populus_euphratica": "Populus euphratica",
        "Populus_simonii": "Populus simonii",
        "Populus_trichocarpa": "Populus trichocarpa",
        "Populus_yunnanensis": "Populus yunnanensis",
    }
    rates = pd.to_numeric(s13["Mapping rate (%)"], errors="coerce")
    for row_index, row in qc.iterrows():
        canonical = species_lookup[str(row["species"])]
        comparison = str(row["comparison"])
        selected = s13[
            s13["Canonical organism"].eq(canonical)
            & s13["Contrast membership"].fillna("").str.contains(comparison, regex=False)
        ]
        selected_rates = rates.loc[selected.index].dropna()
        if selected_rates.empty:
            qc.at[row_index, "mapping_rate_percent"] = pd.NA
            qc.at[row_index, "mapping_rate_note"] = "No retained HISAT2 summary was available for this contrast"
        else:
            qc.at[row_index, "mapping_rate_percent"] = float(selected_rates.median())
            qc.at[row_index, "mapping_rate_note"] = (
                f"Median of archived HISAT2 overall alignment rates from "
                f"{len(selected_rates)}/{len(selected)} libraries"
            )
    return qc


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-root", type=Path, required=True)
    parser.add_argument("--s13", type=Path, required=True)
    parser.add_argument("--contrast-qc", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    run_qc = recover_run_qc(args.archive_root)
    s13 = update_s13(args.s13, run_qc)
    contrast_qc = update_contrast_qc(args.contrast_qc, s13)
    run_qc.to_csv(args.output_dir / "recovered_HISAT2_run_QC.tsv", sep="\t", index=False)
    s13.to_csv(args.s13, sep="\t", index=False, na_rep="NA")
    contrast_qc.to_csv(args.contrast_qc, sep="\t", index=False, na_rep="NA")
    print(
        json.dumps(
            {
                "recovered_run_n": len(run_qc),
                "S13_rows_with_mapping_rate": int(
                    pd.to_numeric(s13["Mapping rate (%)"], errors="coerce").notna().sum()
                ),
                "contrasts_with_mapping_rate": int(
                    pd.to_numeric(contrast_qc["mapping_rate_percent"], errors="coerce").notna().sum()
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
