#!/usr/bin/env python3
"""Build audited TBtools inputs for revised Populus trichocarpa bZIP Figure 4b."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


PAIR_RE = re.compile(r"^\s*(\d+)-\s*(\d+):\s+(\S+)\s+(\S+)\s+(\S+)\s*$")
ALIGN_RE = re.compile(
    r"^## Alignment\s+(\d+):\s+score=(\S+)\s+e_value=(\S+)\s+N=(\d+)\s+"
    r"(\S+)&(\S+)\s+(\S+)\s*$"
)
CHR_RE = re.compile(r"^Chr(?:0?[1-9]|1[0-9])$")

DUPLICATION_COLORS = {
    "WGD/segmental": "#5850C8",
    "dispersed": "#1E90FF",
    "proximal": "#F58220",
    "tandem": "#1FAA59",
    "singleton": "#7F8C8D",
}
CHROMOSOME_COLORS = {
    "Chr01": "#FF8C00",
    "Chr02": "#2E8B57",
    "Chr03": "#DC143C",
    "Chr04": "#8A2BE2",
    "Chr05": "#A0522D",
    "Chr06": "#FF69B4",
    "Chr07": "#708090",
    "Chr08": "#BDB76B",
    "Chr09": "#008B8B",
    "Chr10": "#6B8E23",
    "Chr11": "#FFC125",
    "Chr12": "#48D1CC",
    "Chr13": "#EE6363",
    "Chr14": "#9370DB",
    "Chr15": "#CD853F",
    "Chr16": "#DB7093",
    "Chr17": "#4682B4",
    "Chr18": "#9ACD32",
    "Chr19": "#4169E1",
}
LINK_COLOR = "#DC2D2D"


def rgb_code(hex_color: str) -> str:
    value = hex_color.lstrip("#")
    return f"{int(value[0:2], 16)},{int(value[2:4], 16)},{int(value[4:6], 16)}"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--mcscanx-gff", type=Path, required=True)
    parser.add_argument("--collinearity", type=Path, required=True)
    parser.add_argument("--genome-fai", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_assignments(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    selected = [row for row in rows if row["species"] == "Populus_trichocarpa"]
    if len(selected) != 95:
        raise ValueError(f"Expected 95 P. trichocarpa bZIPs, found {len(selected)}")
    if len({row["geneID"] for row in selected}) != 95:
        raise ValueError("P. trichocarpa gene IDs are not unique")
    return selected


def read_mcscanx_gff(path: Path) -> dict[str, tuple[str, int, int]]:
    records: dict[str, tuple[str, int, int]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            chrom, gene_id, start, end = line.rstrip("\n").split("\t")[:4]
            records[gene_id] = (chrom, int(start), int(end))
    return records


def read_chromosomes(path: Path) -> list[tuple[str, int]]:
    chromosomes: list[tuple[str, int]] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if CHR_RE.match(fields[0]):
                chrom_number = int(fields[0][3:])
                chromosomes.append((f"Chr{chrom_number:02d}", int(fields[1])))
    chromosomes.sort(key=lambda item: int(item[0][3:]))
    if len(chromosomes) != 19:
        raise ValueError(f"Expected 19 chromosomes, found {len(chromosomes)}")
    return chromosomes


def read_collinearity(
    path: Path, selected_ids: set[str]
) -> list[dict[str, str | int]]:
    current: dict[str, str] = {}
    retained: dict[tuple[str, str], dict[str, str | int]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            alignment = ALIGN_RE.match(line)
            if alignment:
                current = {
                    "alignment_id": alignment.group(1),
                    "score": alignment.group(2),
                    "alignment_evalue": alignment.group(3),
                    "block_n": alignment.group(4),
                    "block_chr1": alignment.group(5),
                    "block_chr2": alignment.group(6),
                    "orientation": alignment.group(7),
                }
                continue
            pair = PAIR_RE.match(line)
            if not pair:
                continue
            gene1, gene2 = pair.group(3), pair.group(4)
            if gene1 not in selected_ids or gene2 not in selected_ids:
                continue
            key = tuple(sorted((gene1, gene2)))
            record: dict[str, str | int] = {
                **current,
                "pair_index": int(pair.group(2)),
                "gene1": gene1,
                "gene2": gene2,
                "pair_evalue": pair.group(5),
            }
            retained.setdefault(key, record)
    return sorted(
        retained.values(),
        key=lambda row: (int(str(row["alignment_id"])), int(row["pair_index"])),
    )


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)

    assignments = read_assignments(args.assignments)
    coordinates = read_mcscanx_gff(args.mcscanx_gff)
    chromosomes = read_chromosomes(args.genome_fai)

    missing = sorted(row["geneID"] for row in assignments if row["geneID"] not in coordinates)
    if missing:
        raise ValueError(f"Missing MCScanX coordinates for: {', '.join(missing)}")

    gene_rows: list[dict[str, object]] = []
    for row in assignments:
        chrom, start, end = coordinates[row["geneID"]]
        if not CHR_RE.match(chrom):
            raise ValueError(f"Non-chromosomal retained bZIP: {row['geneID']} on {chrom}")
        short_label = row["pangeneID"].replace("PtbZIP.", "")
        gene_rows.append(
            {
                "Gene_ID": row["geneID"],
                "PtbZIP_ID": row["pangeneID"],
                "Display_label": short_label,
                "Chromosome": chrom,
                "Start": start,
                "End": end,
                "Pan_genome_class": row["class"],
                "Duplication_type": row["type"],
                "Source": row["source"],
                "Color": DUPLICATION_COLORS[row["type"]],
            }
        )
    gene_rows.sort(key=lambda row: (int(str(row["Chromosome"])[3:]), int(row["Start"]), str(row["Gene_ID"])))
    by_gene = {str(row["Gene_ID"]): row for row in gene_rows}

    pairs = read_collinearity(args.collinearity, set(by_gene))
    pair_rows: list[dict[str, object]] = []
    for index, pair in enumerate(pairs, start=1):
        left = by_gene[str(pair["gene1"])]
        right = by_gene[str(pair["gene2"])]
        pair_rows.append(
            {
                "Pair_ID": f"P{index:03d}",
                "Alignment_ID": pair["alignment_id"],
                "Block_score": pair["score"],
                "Block_Evalue": pair["alignment_evalue"],
                "Gene1": pair["gene1"],
                "PtbZIP1": left["PtbZIP_ID"],
                "Chr1": left["Chromosome"],
                "Start1": left["Start"],
                "End1": left["End"],
                "Gene2": pair["gene2"],
                "PtbZIP2": right["PtbZIP_ID"],
                "Chr2": right["Chromosome"],
                "Start2": right["Start"],
                "End2": right["End"],
                "Orientation": pair["orientation"],
                "Pair_Evalue": pair["pair_evalue"],
            }
        )

    # Header-free files for direct TBtools imports.
    with (args.output_dir / "01_TBtools_chromosome_lengths.txt").open("w", encoding="utf-8") as handle:
        for chrom, length in chromosomes:
            handle.write(f"{chrom}\t{length}\n")

    with (args.output_dir / "01b_TBtools_chromosome_lengths_RGB.txt").open("w", encoding="utf-8") as handle:
        for chrom, length in chromosomes:
            handle.write(f"{chrom}\t{length}\t{rgb_code(CHROMOSOME_COLORS[chrom])}\n")

    with (args.output_dir / "02_TBtools_gene_locations_chr_gene_start_end.txt").open("w", encoding="utf-8") as handle:
        for row in gene_rows:
            handle.write(f"{row['Chromosome']}\t{row['Gene_ID']}\t{row['Start']}\t{row['End']}\n")

    with (args.output_dir / "03_TBtools_gene_locations_gene_chr_start_end.txt").open("w", encoding="utf-8") as handle:
        for row in gene_rows:
            handle.write(f"{row['Gene_ID']}\t{row['Chromosome']}\t{row['Start']}\t{row['End']}\n")

    with (args.output_dir / "04_TBtools_collinearity_links_coordinates.txt").open("w", encoding="utf-8") as handle:
        for row in pair_rows:
            handle.write(
                f"{row['Chr1']}\t{row['Start1']}\t{row['End1']}\t"
                f"{row['Chr2']}\t{row['Start2']}\t{row['End2']}\n"
            )

    with (args.output_dir / "04b_TBtools_collinearity_links_RED_RGB.txt").open("w", encoding="utf-8") as handle:
        for row in pair_rows:
            handle.write(
                f"{row['Chr1']}\t{row['Start1']}\t{row['End1']}\t"
                f"{row['Chr2']}\t{row['Start2']}\t{row['End2']}\t{rgb_code(LINK_COLOR)}\n"
            )

    with (args.output_dir / "05_TBtools_collinearity_gene_pairs.txt").open("w", encoding="utf-8") as handle:
        for row in pair_rows:
            handle.write(f"{row['Gene1']}\t{row['Gene2']}\n")

    with (args.output_dir / "06_TBtools_gene_label_track.txt").open("w", encoding="utf-8") as handle:
        for row in gene_rows:
            handle.write(f"{row['Chromosome']}\t{row['Display_label']}\t{row['Start']}\t{row['End']}\n")

    with (args.output_dir / "06b_TBtools_gene_features_duplication_RGB.txt").open("w", encoding="utf-8") as handle:
        for row in gene_rows:
            handle.write(
                f"{row['Chromosome']}\t{row['Display_label']}\t{row['Start']}\t{row['End']}\t"
                f"{rgb_code(str(row['Color']))}\n"
            )

    with (args.output_dir / "07_TBtools_duplication_color_track.txt").open("w", encoding="utf-8") as handle:
        for row in gene_rows:
            handle.write(f"{row['Chromosome']}\t{row['Start']}\t{row['End']}\t{row['Color']}\n")

    # MCScanX-compatible reduced files for TBtools modules that accept GFF + collinearity.
    with (args.output_dir / "08_Ptrichocarpa_95_bZIP_MCScanX.gff").open("w", encoding="utf-8") as handle:
        for row in gene_rows:
            handle.write(f"{row['Chromosome']}\t{row['Gene_ID']}\t{row['Start']}\t{row['End']}\n")

    with (args.output_dir / "09_Ptrichocarpa_strict_bZIP_pairs.collinearity").open("w", encoding="utf-8") as handle:
        handle.write("############ Revised P. trichocarpa bZIP-only collinearity ############\n")
        handle.write(f"# Retained bZIP genes: {len(gene_rows)}\n")
        handle.write(f"# Unique bZIP-bZIP pairs: {len(pair_rows)}\n")
        blocks: dict[str, list[dict[str, object]]] = {}
        for row in pair_rows:
            blocks.setdefault(str(row["Alignment_ID"]), []).append(row)
        for index, rows in enumerate(blocks.values()):
            first = rows[0]
            handle.write(
                f"## Alignment {index}: score={first['Block_score']} e_value={first['Block_Evalue']} "
                f"N={len(rows)} {first['Chr1']}&{first['Chr2']} {first['Orientation']}\n"
            )
            for pair_index, row in enumerate(rows):
                handle.write(
                    f"  {index}- {pair_index:2d}:\t{row['Gene1']}\t{row['Gene2']}\t{row['Pair_Evalue']}\n"
                )

    write_tsv(
        args.output_dir / "10_gene_metadata_with_headers.tsv",
        gene_rows,
        [
            "Gene_ID",
            "PtbZIP_ID",
            "Display_label",
            "Chromosome",
            "Start",
            "End",
            "Pan_genome_class",
            "Duplication_type",
            "Source",
            "Color",
        ],
    )
    write_tsv(
        args.output_dir / "11_collinearity_pair_metadata_with_headers.tsv",
        pair_rows,
        [
            "Pair_ID",
            "Alignment_ID",
            "Block_score",
            "Block_Evalue",
            "Gene1",
            "PtbZIP1",
            "Chr1",
            "Start1",
            "End1",
            "Gene2",
            "PtbZIP2",
            "Chr2",
            "Start2",
            "End2",
            "Orientation",
            "Pair_Evalue",
        ],
    )

    chrom_counts = Counter(str(row["Chromosome"]) for row in gene_rows)
    duplication_counts = Counter(str(row["Duplication_type"]) for row in gene_rows)
    linked_genes = {str(row["Gene1"]) for row in pair_rows} | {str(row["Gene2"]) for row in pair_rows}
    linked_duplication_counts = Counter(
        str(row["Duplication_type"]) for row in gene_rows if str(row["Gene_ID"]) in linked_genes
    )
    unlinked_duplication_counts = Counter(
        str(row["Duplication_type"]) for row in gene_rows if str(row["Gene_ID"]) not in linked_genes
    )
    interchromosomal = sum(row["Chr1"] != row["Chr2"] for row in pair_rows)
    intrachromosomal = len(pair_rows) - interchromosomal
    summary = {
        "retained_bzip_genes": len(gene_rows),
        "assigned_oggs": len({str(row["PtbZIP_ID"]) for row in gene_rows}),
        "chromosomes_with_bzip": sum(count > 0 for count in chrom_counts.values()),
        "chr11_bzip": chrom_counts.get("Chr11", 0),
        "chr12_bzip": chrom_counts.get("Chr12", 0),
        "chromosome_counts": dict(sorted(chrom_counts.items(), key=lambda item: int(item[0][3:]))),
        "duplication_counts": dict(duplication_counts),
        "strict_bzip_bzip_pairs": len(pair_rows),
        "genes_in_strict_network": len(linked_genes),
        "interchromosomal_pairs": interchromosomal,
        "intrachromosomal_pairs": intrachromosomal,
        "unlinked_genes": len(gene_rows) - len(linked_genes),
        "linked_gene_duplication_counts": dict(linked_duplication_counts),
        "unlinked_gene_duplication_counts": dict(unlinked_duplication_counts),
    }
    (args.output_dir / "12_Figure4b_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    readme = f"""Figure 4b revised TBtools input package

Authoritative set
- 95 retained P. trichocarpa bZIP genes: 94 annotated models and 1 audited supplementary model.
- {len(pair_rows)} unique strict bZIP-bZIP MCScanX pairs involving {len(linked_genes)} genes.
- Chr11 and Chr12 contain no retained bZIP gene.

Recommended TBtools inputs (colors embedded)
1. Chromosome/ideogram: 01b_TBtools_chromosome_lengths_RGB.txt
2. Gene labels colored by duplication type: 06b_TBtools_gene_features_duplication_RGB.txt
3. Red collinearity links: 04b_TBtools_collinearity_links_RED_RGB.txt
4. Files 01, 04, and 06 are colorless alternatives.

Alternative formats
- Files 02 and 03 provide both common four-column gene-location orders and use unique full gene IDs as keys.
- File 05 is the strict gene-pair list.
- Files 08 and 09 are reduced MCScanX-compatible GFF and collinearity files for TBtools modules that accept those inputs directly.
- Files 10 and 11 contain headers and full identifiers for scientific auditing; do not use them as header-free TBtools tracks.

Display and interpretation
- Display_label uses the revised pangene suffix (for example, CR001); Gene_ID remains the unique key.
- Duplication colors: WGD/segmental {DUPLICATION_COLORS['WGD/segmental']}; dispersed {DUPLICATION_COLORS['dispersed']}; proximal {DUPLICATION_COLORS['proximal']}; tandem {DUPLICATION_COLORS['tandem']}; singleton {DUPLICATION_COLORS['singleton']}.
- Files 01b, 04b, and 06b embed the chromosome, link, and duplication-type RGB codes recovered from the original Figure 4b SVG.
- The two proteins unassigned by OrthoFinder occur in other species and do not affect this P. trichocarpa panel.

Source files
- Revised assignments: {args.assignments}
- MCScanX GFF: {args.mcscanx_gff}
- MCScanX collinearity: {args.collinearity}
- Genome index: {args.genome_fai}
"""
    (args.output_dir / "README_TBtools_Figure4b.txt").write_text(readme, encoding="utf-8")

    checks = {
        "all_passed": (
            len(gene_rows) == 95
            and len({str(row['Gene_ID']) for row in gene_rows}) == 95
            and len(chromosomes) == 19
            and chrom_counts.get("Chr11", 0) == 0
            and chrom_counts.get("Chr12", 0) == 0
            and len(pair_rows) == 72
            and len(linked_genes) == 74
            and linked_duplication_counts == Counter({"WGD/segmental": 74})
            and all(str(row["Gene1"]) in by_gene and str(row["Gene2"]) in by_gene for row in pair_rows)
        ),
        "summary": summary,
        "source_sha256": {
            "assignments": sha256(args.assignments),
            "mcscanx_gff": sha256(args.mcscanx_gff),
            "collinearity": sha256(args.collinearity),
            "genome_fai": sha256(args.genome_fai),
        },
    }
    (args.output_dir / "QA_REPORT.json").write_text(
        json.dumps(checks, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if not checks["all_passed"]:
        raise RuntimeError("Figure 4b TBtools input QA failed")

    manifest = []
    for output in sorted(args.output_dir.iterdir()):
        if output.is_file():
            manifest.append({"file": output.name, "bytes": output.stat().st_size, "sha256": sha256(output)})
    write_tsv(args.output_dir / "SHA256SUMS.tsv", manifest, ["file", "bytes", "sha256"])

    print(json.dumps(checks, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()





