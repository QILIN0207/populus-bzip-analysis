#!/usr/bin/env python3
"""Parse a MEME Suite MAST XML result and summarize motif occupancy.

The script keeps every input sequence from the stable-ID mapping in the
denominator, including sequences omitted from the XML because their sequence
E-value exceeds MAST's display threshold.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Iterable


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mast-xml", type=Path, required=True)
    parser.add_argument("--mapping-tsv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def pct(numerator: int, denominator: int) -> str:
    return f"{100.0 * numerator / denominator:.6f}" if denominator else ""


def main() -> None:
    args = parse_args()
    mast_xml = args.mast_xml.resolve()
    mapping_tsv = args.mapping_tsv.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    mapping_rows = read_tsv(mapping_tsv)
    if not mapping_rows:
        raise ValueError("Mapping TSV is empty")
    mapping_by_id = {row["upload_id"]: row for row in mapping_rows}
    if len(mapping_by_id) != len(mapping_rows):
        raise ValueError("Duplicate upload_id values in mapping TSV")

    root = ET.parse(mast_xml).getroot()
    settings_node = root.find("settings")
    motif_nodes = root.find("motifs")
    sequences_node = root.find("sequences")
    sequence_db = root.find("sequence_dbs/sequence_db")
    if settings_node is None or motif_nodes is None or sequences_node is None or sequence_db is None:
        raise ValueError("MAST XML is missing required settings, motifs, sequence database, or sequences nodes")

    motifs: list[dict[str, Any]] = []
    for index, node in enumerate(motif_nodes.findall("motif"), start=1):
        motifs.append(
            {
                "motif_number": index,
                "motif_id": node.get("id", ""),
                "motif_alt": node.get("alt", ""),
                "motif_width": int(node.get("length", "0")),
                "discovery_nsites": int(node.get("nsites", "0")),
                "discovery_evalue": node.get("evalue", ""),
            }
        )
    if len(motifs) != 10:
        raise ValueError(f"Expected 10 motifs, found {len(motifs)}")

    hit_rows: list[dict[str, Any]] = []
    sequence_results: dict[str, dict[str, Any]] = {}
    for sequence_node in sequences_node.findall("sequence"):
        upload_id = sequence_node.get("name", "")
        if upload_id not in mapping_by_id:
            raise ValueError(f"MAST sequence is absent from mapping: {upload_id}")
        if upload_id in sequence_results:
            raise ValueError(f"Duplicate MAST sequence: {upload_id}")
        score_node = sequence_node.find("score")
        if score_node is None:
            raise ValueError(f"Missing sequence score for {upload_id}")
        motif_counts: Counter[int] = Counter()
        for segment in sequence_node.findall("seg"):
            for hit in segment.findall("hit"):
                motif_number = int(hit.get("idx", "-1")) + 1
                if not 1 <= motif_number <= len(motifs):
                    raise ValueError(f"Invalid motif index for {upload_id}: {motif_number}")
                hit_pvalue = float(hit.get("pvalue", "nan"))
                motif_counts[motif_number] += 1
                mapping = mapping_by_id[upload_id]
                motif = motifs[motif_number - 1]
                hit_rows.append(
                    {
                        "upload_id": upload_id,
                        "original_id": mapping["original_id"],
                        "species": mapping["species"],
                        "new_orthogroup": mapping["new_orthogroup"],
                        "PtbZIP_ID": mapping["PtbZIP_ID"],
                        "pan_genome_class": mapping["pan_genome_class"],
                        "subfamily": mapping["subfamily"],
                        **motif,
                        "hit_position_1based": int(hit.get("pos", "0")),
                        "hit_pvalue": f"{hit_pvalue:.12g}",
                        "match": hit.get("match", ""),
                        "sequence_combined_pvalue": score_node.get("combined_pvalue", ""),
                        "sequence_evalue": score_node.get("evalue", ""),
                    }
                )
        sequence_results[upload_id] = {
            "sequence_combined_pvalue": score_node.get("combined_pvalue", ""),
            "sequence_evalue": score_node.get("evalue", ""),
            "motif_counts": motif_counts,
        }

    gene_rows: list[dict[str, Any]] = []
    for mapping in mapping_rows:
        upload_id = mapping["upload_id"]
        result = sequence_results.get(upload_id)
        motif_counts = result["motif_counts"] if result else Counter()
        row: dict[str, Any] = dict(mapping)
        row.update(
            {
                "mast_reported": "yes" if result else "no",
                "sequence_combined_pvalue": result["sequence_combined_pvalue"] if result else "",
                "sequence_evalue": result["sequence_evalue"] if result else ">10",
                "total_motif_hits": sum(motif_counts.values()),
                "distinct_motifs_present": len(motif_counts),
                "motifs_present": ",".join(f"motif{number}" for number in sorted(motif_counts)),
            }
        )
        for motif_number in range(1, 11):
            row[f"motif{motif_number}_present"] = 1 if motif_counts[motif_number] else 0
            row[f"motif{motif_number}_hit_count"] = motif_counts[motif_number]
        gene_rows.append(row)

    def summarize(groups: dict[str, list[dict[str, Any]]], group_field: str) -> list[dict[str, Any]]:
        output: list[dict[str, Any]] = []
        for group, rows in groups.items():
            for motif in motifs:
                number = motif["motif_number"]
                genes_with_hit = sum(int(row[f"motif{number}_present"]) for row in rows)
                total_hits = sum(int(row[f"motif{number}_hit_count"]) for row in rows)
                output.append(
                    {
                        group_field: group,
                        **motif,
                        "denominator_genes": len(rows),
                        "genes_with_hit": genes_with_hit,
                        "frequency_percent": pct(genes_with_hit, len(rows)),
                        "total_hits": total_hits,
                    }
                )
        return output

    overall_rows = summarize({"all_primary_input": gene_rows}, "analysis_set")
    assigned_rows = [row for row in gene_rows if row["assignment_status"] == "assigned"]
    subfamily_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    class_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    species_groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in assigned_rows:
        subfamily_groups[row["subfamily"]].append(row)
        class_groups[row["pan_genome_class"]].append(row)
        species_groups[row["species"]].append(row)

    subfamily_order = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "S", "Un"]
    ordered_subfamilies = {key: subfamily_groups[key] for key in subfamily_order}
    class_order = ["core", "softcore", "shell", "cloud"]
    ordered_classes = {key: class_groups[key] for key in class_order}
    ordered_species = {key: species_groups[key] for key in sorted(species_groups)}

    subfamily_summary = summarize(ordered_subfamilies, "subfamily")
    class_summary = summarize(ordered_classes, "pan_genome_class")
    species_summary = summarize(ordered_species, "species")

    group_d = ordered_subfamilies["D"]
    group_d_audit = []
    for row in group_d:
        group_d_audit.append(
            {
                "upload_id": row["upload_id"],
                "original_id": row["original_id"],
                "species": row["species"],
                "new_orthogroup": row["new_orthogroup"],
                "PtbZIP_ID": row["PtbZIP_ID"],
                "pan_genome_class": row["pan_genome_class"],
                "sequence_length": row["sequence_length"],
                "mast_reported": row["mast_reported"],
                "motif3_present": row["motif3_present"],
                "motif3_hit_count": row["motif3_hit_count"],
                "motif3_min_hit_pvalue": min(
                    (
                        float(hit["hit_pvalue"])
                        for hit in hit_rows
                        if hit["upload_id"] == row["upload_id"] and hit["motif_number"] == 3
                    ),
                    default="",
                ),
            }
        )

    omitted_ids = sorted(set(mapping_by_id) - set(sequence_results))
    max_hit_pvalue = max((float(row["hit_pvalue"]) for row in hit_rows), default=0.0)
    motif3_d_count = sum(int(row["motif3_present"]) for row in group_d)
    qa = {
        "all_passed": True,
        "mast_version": root.get("version", ""),
        "mast_release": root.get("release", ""),
        "mast_xml_sha256": sha256(mast_xml),
        "mapping_tsv_sha256": sha256(mapping_tsv),
        "sequence_database": dict(sequence_db.attrib),
        "settings": dict(settings_node.attrib),
        "mapping_rows": len(mapping_rows),
        "assigned_rows": len(assigned_rows),
        "unassigned_rows": len(mapping_rows) - len(assigned_rows),
        "xml_reported_sequences": len(sequence_results),
        "xml_omitted_sequences": len(omitted_ids),
        "omitted_upload_ids": omitted_ids,
        "motifs": len(motifs),
        "total_hits": len(hit_rows),
        "maximum_reported_hit_pvalue": max_hit_pvalue,
        "group_D_denominator": len(group_d),
        "group_D_motif3_genes": motif3_d_count,
        "group_D_motif3_frequency_percent": float(pct(motif3_d_count, len(group_d))),
    }
    expected_input_count = int(sequence_db.get("seq_count", "0"))
    checks = {
        "mapping_matches_database_count": len(mapping_rows) == expected_input_count == 1764,
        "reported_plus_omitted_matches_input": len(sequence_results) + len(omitted_ids) == len(mapping_rows),
        "all_reported_ids_mapped": set(sequence_results).issubset(mapping_by_id),
        "ten_motifs": len(motifs) == 10,
        "hit_pvalues_within_xml_threshold": max_hit_pvalue <= float(settings_node.get("max_hit_pvalue", "nan")),
        "expected_primary_assignment_counts": len(assigned_rows) == 1762 and len(mapping_rows) - len(assigned_rows) == 2,
        "group_D_denominator_nonzero": len(group_d) > 0,
    }
    qa["checks"] = checks
    qa["all_passed"] = all(checks.values())
    if not qa["all_passed"]:
        raise ValueError(f"QA checks failed: {checks}")

    snapshot_path = output_dir / "New1764_primary_MAST_original_result.xml"
    shutil.copy2(mast_xml, snapshot_path)

    mapping_fields = list(mapping_rows[0])
    gene_fields = mapping_fields + [
        "mast_reported",
        "sequence_combined_pvalue",
        "sequence_evalue",
        "total_motif_hits",
        "distinct_motifs_present",
        "motifs_present",
    ] + [field for number in range(1, 11) for field in (f"motif{number}_present", f"motif{number}_hit_count")]
    hit_fields = [
        "upload_id",
        "original_id",
        "species",
        "new_orthogroup",
        "PtbZIP_ID",
        "pan_genome_class",
        "subfamily",
        "motif_number",
        "motif_id",
        "motif_alt",
        "motif_width",
        "discovery_nsites",
        "discovery_evalue",
        "hit_position_1based",
        "hit_pvalue",
        "match",
        "sequence_combined_pvalue",
        "sequence_evalue",
    ]
    summary_fields = [
        "analysis_set",
        "subfamily",
        "pan_genome_class",
        "species",
        "motif_number",
        "motif_id",
        "motif_alt",
        "motif_width",
        "discovery_nsites",
        "discovery_evalue",
        "denominator_genes",
        "genes_with_hit",
        "frequency_percent",
        "total_hits",
    ]
    d_fields = list(group_d_audit[0])

    write_tsv(output_dir / "mast_gene_motif_summary.tsv", gene_rows, gene_fields)
    write_tsv(output_dir / "mast_hits_long.tsv", hit_rows, hit_fields)
    write_tsv(output_dir / "motif_frequency_overall.tsv", overall_rows, summary_fields)
    write_tsv(output_dir / "motif_frequency_by_subfamily.tsv", subfamily_summary, summary_fields)
    write_tsv(output_dir / "motif_frequency_by_pangenome_class.tsv", class_summary, summary_fields)
    write_tsv(output_dir / "motif_frequency_by_species.tsv", species_summary, summary_fields)
    write_tsv(output_dir / "group_D_motif3_audit.tsv", group_d_audit, d_fields)
    (output_dir / "QA_REPORT.json").write_text(json.dumps(qa, indent=2) + "\n", encoding="utf-8")

    print(json.dumps(qa, indent=2))


if __name__ == "__main__":
    main()
