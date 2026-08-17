#!/usr/bin/env python3
"""Rebuild revised Populus bZIP RNA-seq tables after candidate/OGG curation."""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


SPECIES_ALIASES = {
    "Populus_alba_var._pyramidalis": "Populus_alba_var__pyramidalis",
    "Populus_alba_var_pyramidalis": "Populus_alba_var__pyramidalis",
    "Populus_deltoide": "Populus_deltoides",
}

REPLICATED_COMPARISONS = {
    ("Populus_alba_var__pyramidalis", "salt150_vs_control"),
    ("Populus_alba_var__pyramidalis", "salt300_vs_control"),
    ("Populus_deltoides", "salt_vs_control"),
    ("Populus_euphratica", "salt_time12_vs_control_time12"),
    ("Populus_euphratica", "salt_time2_vs_control_time2"),
    ("Populus_simonii", "salt_vs_control"),
    ("Populus_trichocarpa", "prolonged_salt_vs_control"),
    ("Populus_trichocarpa", "short_salt_vs_control"),
}


def normalize_species(value: object) -> str:
    text = str(value).strip()
    return SPECIES_ALIASES.get(text, text)


def normalize_gene_id(value: object) -> str:
    text = str(value).strip()
    if "__" in text:
        text = text.split("__", 1)[1]
    text = re.sub(r"^gene_", "", text)
    text = re.sub(r"\.final$", "", text)
    text = re.sub(r"\.v\d+\.\d+$", "", text)
    text = re.sub(r"\.t\d+$", "", text) if text.startswith("Potri.") else text
    return text


def id_aliases(value: object) -> set[str]:
    base = normalize_gene_id(value)
    aliases = {base}
    aliases.add(re.sub(r"_2dom$", "", base))
    aliases.add(re.sub(r"\.\d+$", "", base))
    if base.startswith("Potri."):
        aliases.add(re.sub(r"\.\d+$", "", base))
    return {x for x in aliases if x}


def truthy(value: object) -> bool:
    return str(value).strip().lower() in {"yes", "true", "1"}


def load_mapping(mapping_path: Path, subfamily_path: Path) -> pd.DataFrame:
    mapping = pd.read_csv(mapping_path, sep="\t")
    mapping["species_norm"] = mapping["species"].map(normalize_species)
    mapping["gene_norm"] = mapping["original_id"].map(normalize_gene_id)
    mapping["pan_genome_class"] = mapping["pan_genome_class"].replace({"softcore": "soft-core"})

    sub = pd.read_csv(subfamily_path, sep="\t")
    sub_cols = {c.lower(): c for c in sub.columns}
    ptb_col = sub_cols.get("ptbzip_id", "PtbZIP_ID")
    sf_col = sub_cols.get("subfamily", "subfamily")
    sf = sub[[ptb_col, sf_col]].drop_duplicates(ptb_col).set_index(ptb_col)[sf_col]
    mapping["subfamily"] = mapping["PtbZIP_ID"].map(sf).fillna(mapping["subfamily"])
    mapping.loc[mapping['PtbZIP_ID'] == 'PtbZIP.CR003', 'subfamily'] = 'S'
    return mapping


def make_alias_lookup(mapping: pd.DataFrame) -> dict[tuple[str, str], list[int]]:
    lookup: dict[tuple[str, str], list[int]] = {}
    for idx, row in mapping.iterrows():
        for alias in id_aliases(row["gene_norm"]):
            lookup.setdefault((row["species_norm"], alias), []).append(idx)
    return lookup


def map_rows_to_revised(
    frame: pd.DataFrame,
    mapping: pd.DataFrame,
    species_col: str,
    gene_cols: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lookup = make_alias_lookup(mapping)
    mapped_rows: list[dict[str, object]] = []
    unmatched_rows: list[dict[str, object]] = []
    for _, source in frame.iterrows():
        species = normalize_species(source[species_col])
        candidates: set[int] = set()
        for col in gene_cols:
            if col not in frame.columns or pd.isna(source[col]):
                continue
            for alias in id_aliases(source[col]):
                candidates.update(lookup.get((species, alias), []))
        if len(candidates) != 1:
            record = source.to_dict()
            record["mapping_candidate_n"] = len(candidates)
            unmatched_rows.append(record)
            continue
        target = mapping.loc[next(iter(candidates))]
        record = source.to_dict()
        record.update(
            {
                "species": species,
                "orthogroup_id": target["new_orthogroup"],
                "pangene_id": target["PtbZIP_ID"],
                "pangenome_class": target["pan_genome_class"],
                "subfamily": target["subfamily"],
                "revised_original_id": target["original_id"],
            }
        )
        mapped_rows.append(record)
    return pd.DataFrame(mapped_rows), pd.DataFrame(unmatched_rows)


def rebuild_s14(
    source_path: Path,
    audit_path: Path,
    mapping: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_excel(source_path, header=1)
    old_by_key: dict[tuple[str, str], pd.Series] = {}
    for _, row in source.iterrows():
        key = (normalize_species(row["species"]), normalize_gene_id(row["expression_gene_id"]))
        old_by_key[key] = row

    audit = pd.read_excel(audit_path)
    retained = audit[audit["Retained in revised bZIP set"].map(truthy)].copy()
    lookup = make_alias_lookup(mapping)
    expression_cols = list(source.columns[6:])

    rows: list[dict[str, object]] = []
    unresolved: list[dict[str, object]] = []
    for _, candidate in retained.iterrows():
        species = normalize_species(candidate["Species"])
        cid = normalize_gene_id(candidate["Candidate ID"])
        map_hits: set[int] = set()
        for alias in id_aliases(cid):
            map_hits.update(lookup.get((species, alias), []))
        if len(map_hits) != 1:
            unresolved.append({"species": species, "candidate_id": cid, "mapping_candidate_n": len(map_hits)})
            continue
        target = mapping.loc[next(iter(map_hits))]
        old = None
        for key, row in old_by_key.items():
            if key[0] == species and bool(id_aliases(key[1]) & id_aliases(cid)):
                old = row
                break
        record: dict[str, object] = {
            "species": species,
            "expression_gene_id": cid,
            "orthogroup_id": target["new_orthogroup"],
            "pangene_id": target["PtbZIP_ID"],
            "subfamily": target["subfamily"],
            "pangenome_class": target["pan_genome_class"],
            "model_classification": candidate["Final classification"],
        }
        for col in expression_cols:
            record[col] = old[col] if old is not None else np.nan
        record.update(
            {
                "samples_evaluated_n": candidate["RNA-seq samples evaluated"],
                "samples_TPM_gt_0_n": candidate["Samples with TPM > 0"],
                "samples_TPM_ge_1_n": candidate["Samples with TPM >= 1"],
                "maximum_TPM": candidate["Maximum TPM"],
                "RNA_seq_support_class": candidate["RNA-seq support class"],
                "public_transcript_evidence": candidate["Public transcript evidence assessment"],
            }
        )
        rows.append(record)
    result = pd.DataFrame(rows)
    order = [
        "species", "expression_gene_id", "orthogroup_id", "pangene_id", "subfamily",
        "pangenome_class", "model_classification", *expression_cols,
        "samples_evaluated_n", "samples_TPM_gt_0_n", "samples_TPM_ge_1_n",
        "maximum_TPM", "RNA_seq_support_class", "public_transcript_evidence",
    ]
    result = result[order].sort_values(["species", "expression_gene_id"]).reset_index(drop=True)
    return result, pd.DataFrame(unresolved)


def rebuild_s15(source_path: Path, mapping: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_excel(source_path, header=1)
    mapped, unmatched = map_rows_to_revised(
        source,
        mapping,
        species_col="species",
        gene_cols=["expression_gene_id", "catalog_gene_id"],
    )
    mapped["replicated_contrast"] = [
        (sp, comp) in REPLICATED_COMPARISONS
        for sp, comp in zip(mapped["species"], mapped["comparison"])
    ]
    mapped["statistically_significant"] = pd.NA
    replicated = mapped["replicated_contrast"]
    mapped.loc[replicated, "statistically_significant"] = (
        mapped.loc[replicated, "log2FoldChange"].abs().ge(1)
        & mapped.loc[replicated, "padj"].le(0.05)
    )
    mapped["descriptive_effect_flag"] = mapped["log2FoldChange"].abs().ge(1)
    mapped["response_direction"] = np.where(
        mapped["log2FoldChange"].isna(),
        pd.NA,
        np.where(mapped["log2FoldChange"] > 0, "up", np.where(mapped["log2FoldChange"] < 0, "down", "no_change")),
    )
    mapped["responsive"] = mapped["statistically_significant"]
    mapped.loc[~replicated, "regulation"] = pd.NA
    out_cols = [
        "species", "comparison", "analysis_type", "catalog_gene_id", "expression_gene_id",
        "orthogroup_id", "pangene_id", "pangenome_class", "subfamily", "log2FoldChange",
        "padj", "replicated_contrast", "statistically_significant", "response_direction",
        "descriptive_effect_flag", "responsive",
    ]
    return mapped[out_cols].reset_index(drop=True), unmatched


def aggregate_contrasts(s15: pd.DataFrame) -> pd.DataFrame:

    rows: list[dict[str, object]] = []
    replicated = s15[s15["replicated_contrast"]].copy()
    for (ogg, ptb, species, comparison), group in replicated.groupby(
        ["orthogroup_id", "pangene_id", "species", "comparison"], dropna=False
    ):
        tested = group[group["log2FoldChange"].notna() & group["padj"].notna()]
        sig = group[group["statistically_significant"].fillna(False).astype(bool)]
        up = sig[sig["log2FoldChange"] > 0]
        down = sig[sig["log2FoldChange"] < 0]
        if len(sig) == 0:
            state = "not_responsive"
        elif len(up) and len(down):
            state = "mixed"
        elif len(up):
            state = "up"
        else:
            state = "down"
        rows.append(
            {
                "orthogroup_id": ogg,
                "pangene_id": ptb,
                "species": species,
                "comparison": comparison,
                "tested_member_n": len(tested),
                "responsive_member_n": len(sig),
                "up_member_n": len(up),
                "down_member_n": len(down),
                "response_state": state,
                "responsive_gene_ids": "; ".join(sorted(sig["catalog_gene_id"].astype(str).tolist())) or pd.NA,
            }
        )
    return pd.DataFrame(rows)


def aggregate_species(contrast_states: pd.DataFrame) -> pd.DataFrame:

    rows: list[dict[str, object]] = []
    for (ogg, ptb, species), group in contrast_states.groupby(
        ["orthogroup_id", "pangene_id", "species"], dropna=False
    ):
        responsive = group[group["response_state"] != "not_responsive"]
        states = set(responsive["response_state"])
        if not states:
            state = "not_responsive"
        elif states == {"up"}:
            state = "up"
        elif states == {"down"}:
            state = "down"
        else:
            state = "mixed"
        rows.append(
            {
                "orthogroup_id": ogg,
                "pangene_id": ptb,
                "species": species,
                "replicated_contrast_n": len(group),
                "responsive_contrast_n": len(responsive),
                "responsive_member_contrast_n": int(responsive["responsive_member_n"].sum()),
                "species_response_state": state,
                "responsive_comparisons": "; ".join(responsive["comparison"].tolist()) or pd.NA,
            }
        )
    return pd.DataFrame(rows)


def aggregate_cross_species(species_states: pd.DataFrame) -> pd.DataFrame:

    rows: list[dict[str, object]] = []
    for (ogg, ptb), group in species_states.groupby(["orthogroup_id", "pangene_id"], dropna=False):
        responsive = group[group["species_response_state"] != "not_responsive"]
        up_n = int((responsive["species_response_state"] == "up").sum())
        down_n = int((responsive["species_response_state"] == "down").sum())
        mixed_n = int((responsive["species_response_state"] == "mixed").sum())
        concordant = (up_n >= 2 and down_n == 0 and mixed_n == 0) or (down_n >= 2 and up_n == 0 and mixed_n == 0)
        if len(responsive) == 0:
            pattern = "not_responsive"
        elif concordant:
            pattern = "directionally_concordant_multiple_species"
        elif len(responsive) == 1:
            pattern = "single_species_only"
        else:
            pattern = "multiple_species_nonconcordant"
        rows.append(
            {
                "orthogroup_id": ogg,
                "pangene_id": ptb,
                "species_tested_n": len(group),
                "responsive_species_n": len(responsive),
                "up_species_n": up_n,
                "down_species_n": down_n,
                "mixed_species_n": mixed_n,
                "cross_species_pattern": pattern,
                "responsive_species": "; ".join(
                    f"{r.species}:{r.species_response_state}" for r in responsive.itertuples()
                ) or pd.NA,
            }
        )
    return pd.DataFrame(rows)


def build_contrast_qc(s15: pd.DataFrame, runs_path: Path) -> pd.DataFrame:
    runs = pd.read_csv(runs_path, sep="\t")
    study_for_species = {
        "Populus_alba_var__pyramidalis": "PRJNA393495",
        "Populus_deltoides": "PRJNA952677",
        "Populus_euphratica": "PRJEB37975",
        "Populus_simonii": "PRJNA359403",
        "Populus_trichocarpa": "PRJEB19784",
        "Populus_yunnanensis": "PRJNA1222559",
    }
    euphratica_runs = {
        "salt_time12_vs_control_time12": {
            "ERR4059418", "ERR4059421", "ERR4059424", "ERR4059427", "ERR4059430",
            "ERR4059448", "ERR4059454", "ERR4059460", "ERR4059472", "ERR4059475",
        },
        "salt_time2_vs_control_time2": {
            "ERR4059433", "ERR4059436", "ERR4059439", "ERR4059442", "ERR4059445",
            "ERR4059451", "ERR4059457", "ERR4059463", "ERR4059466", "ERR4059469",
        },
    }
    rows: list[dict[str, object]] = []
    for (species, comparison), group in s15.groupby(["species", "comparison"]):
        study = study_for_species[species]
        study_runs = runs[runs["study_accession"] == study]
        if species == "Populus_alba_var__pyramidalis":
            if comparison.startswith("salt150"):
                selected = study_runs[study_runs["experiment_title"].str.contains("under 150 mM|under 0 mM", regex=True)]
            elif comparison.startswith("salt300"):
                selected = study_runs[study_runs["experiment_title"].str.contains("under 300 mM|under 0 mM", regex=True)]
            else:
                selected = study_runs
        elif species == "Populus_trichocarpa":
            if comparison.startswith("short"):
                selected = study_runs[study_runs["sample_title"].str.contains("control|short-term", case=False, regex=True)]
            else:
                selected = study_runs[study_runs["sample_title"].str.contains("control|prolonged", case=False, regex=True)]
        elif species == "Populus_euphratica":
            # Run sets follow the official E-MTAB-8988 SDRF assignments.
            selected = study_runs[study_runs["run_accession"].isin(euphratica_runs[comparison])]
        elif species == "Populus_yunnanensis":
            label = "Treat1|Control" if comparison.startswith("treat1") else "Treat4|Control"
            selected = study_runs[study_runs["experiment_title"].str.contains(label, case=False, regex=True)]
        else:
            selected = study_runs
        replicated = (species, comparison) in REPLICATED_COMPARISONS
        rows.append(
            {
                "species": species,
                "comparison": comparison,
                "study_accession": study,
                "analysis_type": "DESeq2" if replicated else "descriptive_TPM_fold_change",
                "biological_replicates_per_group": 3 if replicated and species != "Populus_euphratica" else (5 if replicated else 1),
                "library_n": len(selected),
                "median_raw_reads": float(selected["read_count"].median()),
                "median_raw_bases": float(selected["base_count"].median()),
                "mapping_rate_percent": pd.NA,
                "mapping_rate_note": "Not recoverable from retained historical alignment logs",
                "catalog_bzip_rows_n": len(group),
                "testable_bzip_n": int((group["log2FoldChange"].notna() & group["padj"].notna()).sum()) if replicated else int(group["log2FoldChange"].notna().sum()),
                "significant_response_n": int(group["statistically_significant"].fillna(False).astype(bool).sum()) if replicated else pd.NA,
                "descriptive_abs_log2FC_ge_1_n": int(group["descriptive_effect_flag"].fillna(False).astype(bool).sum()),
                "control_note": "Shared untreated control used for both short and prolonged contrasts" if species == "Populus_trichocarpa" else pd.NA,
            }
        )
    return pd.DataFrame(rows)


def write_outputs(args: argparse.Namespace) -> None:
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    mapping = load_mapping(args.mapping, args.subfamily_mapping)
    s14, s14_unresolved = rebuild_s14(args.s14, args.candidate_audit, mapping)
    s15, s15_unmatched = rebuild_s15(args.s15, mapping)
    contrast = aggregate_contrasts(s15)
    species = aggregate_species(contrast)
    cross = aggregate_cross_species(species)
    qc = build_contrast_qc(s15, args.selected_runs)

    outputs = {
        "S14_revised.tsv": s14,
        "S14_unresolved.tsv": s14_unresolved,
        "S15_revised.tsv": s15,
        "S15_removed_unmatched.tsv": s15_unmatched,
        "S16_contrast_level.tsv": contrast,
        "S16_species_level.tsv": species,
        "S16_cross_species_summary.tsv": cross,
        "contrast_QC.tsv": qc,
    }
    for name, frame in outputs.items():
        frame.to_csv(output_dir / name, sep="\t", index=False, na_rep="NA")

    replicated = s15[s15["replicated_contrast"]]
    significant = replicated[replicated["statistically_significant"].fillna(False).astype(bool)]
    summary = {
        "s14_retained_models": len(s14),
        "s14_evaluable_models": int(s14["samples_evaluated_n"].notna().sum()),
        "s14_tpm_ge_1_models": int((pd.to_numeric(s14["samples_TPM_ge_1_n"], errors="coerce") > 0).sum()),
        "s14_unresolved_models": len(s14_unresolved),
        "s15_revised_rows": len(s15),
        "s15_removed_rows": len(s15_unmatched),
        "replicated_significant_events": len(significant),
        "replicated_significant_unique_genes": int(significant["catalog_gene_id"].nunique()),
        "replicated_significant_oggs": int(significant["orthogroup_id"].nunique()),
        "up_events": int((significant["log2FoldChange"] > 0).sum()),
        "down_events": int((significant["log2FoldChange"] < 0).sum()),
        "cross_species_concordant_oggs": int((cross["cross_species_pattern"] == "directionally_concordant_multiple_species").sum()),
        "cross_species_nonconcordant_multiple_oggs": int((cross["cross_species_pattern"] == "multiple_species_nonconcordant").sum()),
        "single_species_responsive_oggs": int((cross["cross_species_pattern"] == "single_species_only").sum()),
    }
    (output_dir / "RNAseq_revision_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--s14", type=Path, required=True)
    parser.add_argument("--s15", type=Path, required=True)
    parser.add_argument("--candidate-audit", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--subfamily-mapping", type=Path, required=True)
    parser.add_argument("--selected-runs", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    write_outputs(parse_args())




