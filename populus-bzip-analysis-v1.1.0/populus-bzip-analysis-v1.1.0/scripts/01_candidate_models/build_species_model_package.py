#!/usr/bin/env python3
"""Build clean, standard, per-species GFF3/CDS/protein deliverables.

The output GFF3 files contain only the standard gene -> mRNA -> exon/CDS
hierarchy. Audit metadata is kept in a separate ID-mapping table rather than in
GFF3 attributes. Sequence identifiers are cleaned consistently across GFF3,
CDS FASTA, and protein FASTA.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

SOURCE = "bZIP_reannotation"


@dataclass
class Feature:
    seqid: str
    feature_type: str
    start: int
    end: int
    strand: str
    phase: str
    attributes: dict[str, str]
    model_id: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-gff", required=True, type=Path)
    parser.add_argument("--cds-fasta", required=True, type=Path)
    parser.add_argument("--protein-fasta", required=True, type=Path)
    parser.add_argument("--audit-tsv", required=True, type=Path)
    parser.add_argument("--species-list", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def read_gff(path: Path) -> list[Feature]:
    records: list[Feature] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise RuntimeError(f"Expected 9 GFF3 columns: {line}")
            records.append(
                Feature(
                    seqid=fields[0],
                    feature_type=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    strand=fields[6],
                    phase=fields[7],
                    attributes=parse_attributes(fields[8]),
                )
            )
    return records


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    identifier: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    records[identifier] = "".join(chunks).upper()
                identifier = line[1:].split()[0]
                chunks = []
            else:
                if identifier is None:
                    raise RuntimeError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if identifier is not None:
        records[identifier] = "".join(chunks).upper()
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, rows: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in rows:
            handle.write(f">{identifier}\n")
            for index in range(0, len(sequence), 60):
                handle.write(sequence[index : index + 60] + "\n")


def translate(sequence: str) -> str:
    return "".join(
        CODON_TABLE.get(sequence[index : index + 3], "X")
        for index in range(0, len(sequence) - 2, 3)
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def attrs_text(attributes: dict[str, str]) -> str:
    return ";".join(f"{key}={value}" for key, value in attributes.items()) + ";"


def gff_line(
    seqid: str,
    feature_type: str,
    start: int,
    end: int,
    strand: str,
    phase: str,
    attributes: dict[str, str],
) -> str:
    return "\t".join(
        [
            seqid,
            SOURCE,
            feature_type,
            str(start),
            str(end),
            ".",
            strand,
            phase,
            attrs_text(attributes),
        ]
    )


def assign_models(features: list[Feature]) -> tuple[dict[str, Feature], dict[str, Feature]]:
    gene_id_to_model: dict[str, str] = {}
    mrna_id_to_model: dict[str, str] = {}
    genes: dict[str, Feature] = {}
    mrnas: dict[str, Feature] = {}
    for feature in features:
        if feature.feature_type == "gene":
            model_id = feature.attributes.get(
                "Name", feature.attributes["ID"].removesuffix(".gene")
            )
            feature.model_id = model_id
            gene_id_to_model[feature.attributes["ID"]] = model_id
            genes[model_id] = feature
    for feature in features:
        if feature.feature_type in {"mRNA", "transcript"}:
            parent = feature.attributes.get("Parent", "")
            if parent not in gene_id_to_model:
                raise RuntimeError(f"mRNA parent absent: {feature.attributes.get('ID')}")
            model_id = gene_id_to_model[parent]
            feature.model_id = model_id
            mrna_id_to_model[feature.attributes["ID"]] = model_id
            mrnas[model_id] = feature
    for feature in features:
        if feature.feature_type in {"exon", "CDS"}:
            parent = feature.attributes.get("Parent", "")
            if parent not in mrna_id_to_model:
                raise RuntimeError(f"Child parent absent: {feature.attributes.get('ID')}")
            feature.model_id = mrna_id_to_model[parent]
    if set(genes) != set(mrnas):
        raise RuntimeError("Gene and mRNA model sets differ")
    return genes, mrnas


def main() -> None:
    args = parse_args()
    for path in [
        args.input_gff,
        args.cds_fasta,
        args.protein_fasta,
        args.audit_tsv,
        args.species_list,
    ]:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output_dir}")

    features = read_gff(args.input_gff)
    genes, mrnas = assign_models(features)
    cds_sequences = read_fasta(args.cds_fasta)
    protein_sequences = read_fasta(args.protein_fasta)
    audit_rows = read_tsv(args.audit_tsv)
    retained_rows = [row for row in audit_rows if row["retained_in_revised_bZIP_set"] == "Yes"]
    audit_by_model = {row["deliverable_model_id"]: row for row in retained_rows}
    retained_ids = set(audit_by_model)
    if not (
        set(genes)
        == set(mrnas)
        == set(cds_sequences)
        == set(protein_sequences)
        == retained_ids
    ):
        raise RuntimeError("GFF3, FASTA, and retained audit ID sets differ")

    species_rows = read_tsv(args.species_list)
    all_species = [row["species"] for row in species_rows]
    if len(all_species) != 19 or len(set(all_species)) != 19:
        raise RuntimeError("Expected 19 unique species in species list")

    exons_by_model: dict[str, list[Feature]] = defaultdict(list)
    cds_by_model: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        if feature.feature_type == "exon":
            exons_by_model[feature.model_id].append(feature)
        elif feature.feature_type == "CDS":
            cds_by_model[feature.model_id].append(feature)

    clean_by_old = {
        old_id: old_id.removesuffix(".final") for old_id in sorted(retained_ids)
    }
    if len(set(clean_by_old.values())) != 25:
        raise RuntimeError("Removing .final created duplicate model IDs")
    if any("final" in identifier.lower() for identifier in clean_by_old.values()):
        raise RuntimeError("A cleaned model ID still contains final")

    lines_by_species: dict[str, list[str]] = defaultdict(lambda: ["##gff-version 3"])
    cds_by_species: dict[str, list[tuple[str, str]]] = defaultdict(list)
    protein_by_species: dict[str, list[tuple[str, str]]] = defaultdict(list)
    mapping_rows: list[dict[str, object]] = []
    output_ids: list[str] = []
    feature_counts = Counter()

    for old_id in sorted(retained_ids, key=lambda value: (audit_by_model[value]["species"], value)):
        audit = audit_by_model[old_id]
        species = audit["species"]
        if species not in all_species:
            raise RuntimeError(f"Retained model species absent from 19-species list: {species}")
        gene = genes[old_id]
        mrna = mrnas[old_id]
        gene_id = clean_by_old[old_id]
        transcript_id = f"{gene_id}.t1"
        output_ids.append(transcript_id)

        if not (
            gene.seqid == mrna.seqid
            and gene.strand == mrna.strand
            and gene.start == mrna.start
            and gene.end == mrna.end
        ):
            raise RuntimeError(f"Gene/mRNA span mismatch: {old_id}")

        model_lines = [
            gff_line(
                gene.seqid,
                "gene",
                gene.start,
                gene.end,
                gene.strand,
                ".",
                {"ID": gene_id, "Name": gene_id},
            ),
            gff_line(
                mrna.seqid,
                "mRNA",
                mrna.start,
                mrna.end,
                mrna.strand,
                ".",
                {"ID": transcript_id, "Parent": gene_id, "Name": transcript_id},
            ),
        ]
        feature_counts.update(["gene", "mRNA"])

        reverse = gene.strand == "-"
        exons = sorted(exons_by_model[old_id], key=lambda item: item.start, reverse=reverse)
        cds_records = sorted(cds_by_model[old_id], key=lambda item: item.start, reverse=reverse)
        if len(exons) != len(cds_records):
            raise RuntimeError(f"Expected one exon per CDS for {old_id}")

        cumulative = 0
        used_cds: set[int] = set()
        for exon_index, exon in enumerate(exons, start=1):
            exon_id = f"{transcript_id}.exon{exon_index}"
            model_lines.append(
                gff_line(
                    exon.seqid,
                    "exon",
                    exon.start,
                    exon.end,
                    exon.strand,
                    ".",
                    {"ID": exon_id, "Parent": transcript_id},
                )
            )
            feature_counts.update(["exon"])
            contained = [
                cds
                for cds in cds_records
                if id(cds) not in used_cds
                and exon.start <= cds.start
                and cds.end <= exon.end
            ]
            if len(contained) != 1:
                raise RuntimeError(
                    f"Expected exactly one CDS in exon {exon_index} for {old_id}"
                )
            cds = contained[0]
            used_cds.add(id(cds))
            phase = (3 - cumulative % 3) % 3
            cds_id = f"{transcript_id}.cds{exon_index}"
            model_lines.append(
                gff_line(
                    cds.seqid,
                    "CDS",
                    cds.start,
                    cds.end,
                    cds.strand,
                    str(phase),
                    {"ID": cds_id, "Parent": transcript_id},
                )
            )
            feature_counts.update(["CDS"])
            cumulative += cds.end - cds.start + 1
        if len(used_cds) != len(cds_records):
            raise RuntimeError(f"Unmatched CDS records remain for {old_id}")
        if cumulative != len(cds_sequences[old_id]):
            raise RuntimeError(f"GFF3 CDS length differs from FASTA for {old_id}")
        translated = translate(cds_sequences[old_id])
        if not translated.endswith("*") or "*" in translated[:-1]:
            raise RuntimeError(f"Invalid complete ORF for {old_id}")
        if translated[:-1] != protein_sequences[old_id].rstrip("*"):
            raise RuntimeError(f"CDS/protein mismatch for {old_id}")

        lines_by_species[species].extend(model_lines)
        cds_by_species[species].append((transcript_id, cds_sequences[old_id]))
        protein_by_species[species].append(
            (transcript_id, protein_sequences[old_id].rstrip("*"))
        )
        mapping_rows.append(
            {
                "species": species,
                "source_model_id": old_id,
                "standard_gene_id": gene_id,
                "standard_transcript_id": transcript_id,
                "candidate_gene_id": audit["candidate_gene_id"],
                "counting_category": audit["final_counting_category"],
                "independent_new_gene": audit["independent_new_gene"],
            }
        )

    if len(output_ids) != 25 or len(set(output_ids)) != 25:
        raise RuntimeError("Expected 25 unique standardized transcript IDs")
    if feature_counts != Counter({"gene": 25, "mRNA": 25, "exon": 108, "CDS": 108}):
        raise RuntimeError(f"Unexpected output feature counts: {feature_counts}")

    args.output_dir.mkdir(parents=True)
    by_species_root = args.output_dir / "by_species"
    by_species_root.mkdir()
    for species in all_species:
        if species not in lines_by_species:
            continue
        species_dir = by_species_root / species
        species_dir.mkdir()
        stem = f"{species}.supplemented_bZIP_models"
        (species_dir / f"{stem}.gff3").write_text(
            "\n".join(lines_by_species[species]) + "\n", encoding="utf-8"
        )
        write_fasta(species_dir / f"{stem}.cds.fasta", cds_by_species[species])
        write_fasta(
            species_dir / f"{stem}.protein.fasta", protein_by_species[species]
        )

    write_tsv(
        args.output_dir / "model_id_mapping.tsv",
        mapping_rows,
        [
            "species",
            "source_model_id",
            "standard_gene_id",
            "standard_transcript_id",
            "candidate_gene_id",
            "counting_category",
            "independent_new_gene",
        ],
    )

    counts_by_species = Counter(row["species"] for row in mapping_rows)
    independent_by_species = Counter(
        row["species"] for row in mapping_rows if row["independent_new_gene"] == "Yes"
    )
    correction_by_species = Counter(
        row["species"] for row in mapping_rows if row["independent_new_gene"] == "No"
    )
    species_count_rows = [
        {
            "species": species,
            "retained_models": counts_by_species[species],
            "independent_new_loci": independent_by_species[species],
            "annotation_corrections": correction_by_species[species],
            "deliverable_directory": (
                f"by_species/{species}" if counts_by_species[species] else "NA"
            ),
        }
        for species in all_species
    ]
    write_tsv(
        args.output_dir / "species_model_counts.tsv",
        species_count_rows,
        [
            "species",
            "retained_models",
            "independent_new_loci",
            "annotation_corrections",
            "deliverable_directory",
        ],
    )

    gff_text = "\n".join(
        line
        for species in all_species
        for line in lines_by_species.get(species, [])
    )
    forbidden_tokens = [
        ".final",
        "Final52_audit",
        "species=",
        "source_candidate=",
        "orthogroup=",
        "counting_category=",
    ]
    present_forbidden = [token for token in forbidden_tokens if token in gff_text]
    if present_forbidden:
        raise RuntimeError(f"Forbidden GFF3 content remains: {present_forbidden}")

    allowed = {
        "gene": {"ID", "Name"},
        "mRNA": {"ID", "Parent", "Name"},
        "exon": {"ID", "Parent"},
        "CDS": {"ID", "Parent"},
    }
    for species, lines in lines_by_species.items():
        for line in lines[1:]:
            fields = line.split("\t")
            keys = set(parse_attributes(fields[8]))
            if keys != allowed[fields[2]]:
                raise RuntimeError(
                    f"Unexpected attributes in {species} {fields[2]}: {keys}"
                )

    report = {
        "retained_models": 25,
        "independent_new_loci": 21,
        "annotation_corrections": 4,
        "species_with_models": sum(value > 0 for value in counts_by_species.values()),
        "species_without_models": sum(counts_by_species[species] == 0 for species in all_species),
        "feature_counts": dict(feature_counts),
        "gff3_source_column": SOURCE,
        "gff3_attribute_policy": {
            "gene": ["ID", "Name"],
            "mRNA": ["ID", "Parent", "Name"],
            "exon": ["ID", "Parent"],
            "CDS": ["ID", "Parent"],
        },
        "forbidden_tokens_absent": True,
        "all_internal_checks_passed": True,
    }
    (args.output_dir / "STANDARD_GFF3_BUILD_REPORT.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    files = sorted(
        path
        for path in args.output_dir.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.tsv"
    )
    with (args.output_dir / "SHA256SUMS.tsv").open(
        "w", encoding="utf-8", newline="\n"
    ) as handle:
        handle.write("sha256\trelative_path\n")
        for path in files:
            handle.write(f"{sha256(path)}\t{path.relative_to(args.output_dir)}\n")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
