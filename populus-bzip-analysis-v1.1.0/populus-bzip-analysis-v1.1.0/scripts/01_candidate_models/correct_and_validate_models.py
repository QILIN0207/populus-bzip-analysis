#!/usr/bin/env python3
"""Create a corrected, reusable GFF3 package for retained Populus bZIP models.

The script preserves the source package, recalculates CDS phase from the final
complete ORF, adds exon features when an mRNA has CDS features but no exon
features, optionally replaces species-specific sequence aliases, and validates
the corrected hierarchy against the retained CDS/protein FASTA files.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
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


@dataclass
class Feature:
    seqid: str
    source: str
    feature_type: str
    start: int
    end: int
    score: str
    strand: str
    phase: str
    attributes: dict[str, str]
    model_id: str = ""

    def to_line(self) -> str:
        attrs = ";".join(f"{key}={value}" for key, value in self.attributes.items())
        return "\t".join(
            [
                self.seqid,
                self.source,
                self.feature_type,
                str(self.start),
                str(self.end),
                self.score,
                self.strand,
                self.phase,
                attrs,
            ]
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-gff", required=True, type=Path)
    parser.add_argument("--cds-fasta", required=True, type=Path)
    parser.add_argument("--protein-fasta", required=True, type=Path)
    parser.add_argument("--audit-tsv", required=True, type=Path)
    parser.add_argument("--seqid-map", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    return parser.parse_args()


def parse_attributes(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def read_gff(path: Path) -> list[Feature]:
    features: list[Feature] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise RuntimeError(f"Expected 9 GFF3 columns: {line}")
            features.append(
                Feature(
                    seqid=fields[0],
                    source=fields[1],
                    feature_type=fields[2],
                    start=int(fields[3]),
                    end=int(fields[4]),
                    score=fields[5],
                    strand=fields[6],
                    phase=fields[7],
                    attributes=parse_attributes(fields[8]),
                )
            )
    return features


def read_fasta(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    sequences: dict[str, str] = {}
    headers: dict[str, str] = {}
    identifier: str | None = None
    header = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    sequences[identifier] = "".join(chunks).upper()
                    headers[identifier] = header
                header = line[1:]
                identifier = header.split()[0]
                chunks = []
            else:
                if identifier is None:
                    raise RuntimeError(f"Sequence before FASTA header in {path}")
                chunks.append(line)
    if identifier is not None:
        sequences[identifier] = "".join(chunks).upper()
        headers[identifier] = header
    return sequences, headers


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, identifiers: list[str], sequences: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier in identifiers:
            handle.write(f">{identifier}\n")
            sequence = sequences[identifier]
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


def assign_models(features: list[Feature]) -> tuple[dict[str, Feature], dict[str, Feature]]:
    gene_id_to_model: dict[str, str] = {}
    genes: dict[str, Feature] = {}
    mrna_id_to_model: dict[str, str] = {}
    mrnas: dict[str, Feature] = {}

    for feature in features:
        if feature.feature_type == "gene":
            gene_id = feature.attributes["ID"]
            model_id = feature.attributes.get("Name", gene_id.removesuffix(".gene"))
            feature.model_id = model_id
            gene_id_to_model[gene_id] = model_id
            genes[model_id] = feature

    for feature in features:
        if feature.feature_type in {"mRNA", "transcript"}:
            parent = feature.attributes.get("Parent", "")
            if parent not in gene_id_to_model:
                raise RuntimeError(f"mRNA parent absent: {feature.attributes.get('ID')}")
            model_id = gene_id_to_model[parent]
            feature.model_id = model_id
            mrna_id = feature.attributes["ID"]
            mrna_id_to_model[mrna_id] = model_id
            mrnas[model_id] = feature

    for feature in features:
        if feature.feature_type in {"CDS", "exon"}:
            parent = feature.attributes.get("Parent", "")
            if parent not in mrna_id_to_model:
                raise RuntimeError(f"Child parent absent: {feature.attributes.get('ID')}")
            feature.model_id = mrna_id_to_model[parent]

    if set(genes) != set(mrnas):
        raise RuntimeError("Gene and mRNA model sets differ")
    return genes, mrnas


def main() -> None:
    args = parse_args()
    inputs = [args.input_gff, args.cds_fasta, args.protein_fasta, args.audit_tsv, args.seqid_map]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    if args.output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {args.output_dir}")

    features = read_gff(args.input_gff)
    genes, mrnas = assign_models(features)
    cds_sequences, _ = read_fasta(args.cds_fasta)
    protein_sequences, _ = read_fasta(args.protein_fasta)
    audit_rows = read_tsv(args.audit_tsv)
    retained_rows = [row for row in audit_rows if row["retained_in_revised_bZIP_set"] == "Yes"]
    retained_ids = {row["deliverable_model_id"] for row in retained_rows}

    model_sets = [set(genes), set(mrnas), set(cds_sequences), set(protein_sequences), retained_ids]
    if any(current != retained_ids for current in model_sets):
        raise RuntimeError("GFF3, FASTA, and retained-audit model ID sets differ")
    if (len(retained_rows), sum(row["independent_new_gene"] == "Yes" for row in retained_rows)) != (25, 21):
        raise RuntimeError("Expected 25 retained models including 21 independent loci")

    seqid_rows = read_tsv(args.seqid_map)
    seqid_map = {
        (row["species"], row["old_seqid"]): row["new_seqid"] for row in seqid_rows
    }
    species_by_model = {
        model_id: gene.attributes["species"] for model_id, gene in genes.items()
    }

    seqid_changes: list[dict[str, object]] = []
    for feature in features:
        if not feature.model_id:
            continue
        species = species_by_model[feature.model_id]
        mapped = seqid_map.get((species, feature.seqid), feature.seqid)
        if mapped != feature.seqid:
            seqid_changes.append(
                {
                    "model_id": feature.model_id,
                    "feature_id": feature.attributes.get("ID", "NA"),
                    "feature_type": feature.feature_type,
                    "old_seqid": feature.seqid,
                    "new_seqid": mapped,
                }
            )
            feature.seqid = mapped

    cds_by_mrna: dict[str, list[Feature]] = defaultdict(list)
    exons_by_mrna: dict[str, list[Feature]] = defaultdict(list)
    for feature in features:
        if feature.feature_type == "CDS":
            cds_by_mrna[feature.attributes["Parent"]].append(feature)
        elif feature.feature_type == "exon":
            exons_by_mrna[feature.attributes["Parent"]].append(feature)

    phase_changes: list[dict[str, object]] = []
    for mrna_id, records in cds_by_mrna.items():
        reverse = records[0].strand == "-"
        ordered = sorted(records, key=lambda item: item.start, reverse=reverse)
        cumulative = 0
        for index, feature in enumerate(ordered, start=1):
            expected = (3 - cumulative % 3) % 3
            old_phase = feature.phase
            if old_phase != str(expected):
                phase_changes.append(
                    {
                        "model_id": feature.model_id,
                        "mRNA_id": mrna_id,
                        "CDS_index_transcript_order": index,
                        "start": feature.start,
                        "end": feature.end,
                        "strand": feature.strand,
                        "old_phase": old_phase,
                        "new_phase": expected,
                    }
                )
            feature.phase = str(expected)
            cumulative += feature.end - feature.start + 1

    existing_ids = {
        feature.attributes["ID"] for feature in features if "ID" in feature.attributes
    }
    added_exons: list[Feature] = []
    added_exon_rows: list[dict[str, object]] = []
    first_cds_by_mrna: dict[str, Feature] = {}
    for mrna_id, cds_records in cds_by_mrna.items():
        first_cds_by_mrna[mrna_id] = min(cds_records, key=lambda item: features.index(item))
        if exons_by_mrna.get(mrna_id):
            continue
        for index, cds_feature in enumerate(sorted(cds_records, key=lambda item: item.start), start=1):
            exon_id = f"{mrna_id}.exon{index}"
            if exon_id in existing_ids:
                raise RuntimeError(f"Generated exon ID already exists: {exon_id}")
            existing_ids.add(exon_id)
            exon = Feature(
                seqid=cds_feature.seqid,
                source=cds_feature.source,
                feature_type="exon",
                start=cds_feature.start,
                end=cds_feature.end,
                score=cds_feature.score,
                strand=cds_feature.strand,
                phase=".",
                attributes={"ID": exon_id, "Parent": mrna_id},
                model_id=cds_feature.model_id,
            )
            added_exons.append(exon)
            added_exon_rows.append(
                {
                    "model_id": exon.model_id,
                    "mRNA_id": mrna_id,
                    "exon_id": exon_id,
                    "seqid": exon.seqid,
                    "start": exon.start,
                    "end": exon.end,
                    "strand": exon.strand,
                }
            )

    added_by_first_cds: dict[int, list[Feature]] = defaultdict(list)
    for exon in added_exons:
        mrna_id = exon.attributes["Parent"]
        added_by_first_cds[id(first_cds_by_mrna[mrna_id])].append(exon)

    corrected: list[Feature] = []
    for feature in features:
        if id(feature) in added_by_first_cds:
            corrected.extend(sorted(added_by_first_cds[id(feature)], key=lambda item: item.start))
        corrected.append(feature)

    all_feature_ids = [
        feature.attributes["ID"] for feature in corrected if "ID" in feature.attributes
    ]
    duplicates = [identifier for identifier, count in Counter(all_feature_ids).items() if count > 1]
    if duplicates:
        raise RuntimeError(f"Duplicate GFF3 feature IDs: {duplicates}")

    corrected_cds_by_model: dict[str, list[Feature]] = defaultdict(list)
    corrected_exons_by_model: dict[str, list[Feature]] = defaultdict(list)
    for feature in corrected:
        if feature.feature_type == "CDS":
            corrected_cds_by_model[feature.model_id].append(feature)
        elif feature.feature_type == "exon":
            corrected_exons_by_model[feature.model_id].append(feature)

    failures: list[str] = []
    for model_id in sorted(retained_ids):
        cds_records = corrected_cds_by_model[model_id]
        cds_length = sum(feature.end - feature.start + 1 for feature in cds_records)
        if cds_length != len(cds_sequences[model_id]):
            failures.append(f"{model_id}: GFF3 CDS length differs from FASTA")
        translated = translate(cds_sequences[model_id])
        if not translated.endswith("*") or "*" in translated[:-1]:
            failures.append(f"{model_id}: CDS does not encode one terminal stop")
        if translated[:-1] != protein_sequences[model_id].rstrip("*"):
            failures.append(f"{model_id}: CDS translation differs from protein")
        if not corrected_exons_by_model[model_id]:
            failures.append(f"{model_id}: no exon features")
        for cds_feature in cds_records:
            if not any(
                exon.start <= cds_feature.start and cds_feature.end <= exon.end
                for exon in corrected_exons_by_model[model_id]
            ):
                failures.append(f"{model_id}: CDS is not covered by an exon")
        reverse = cds_records[0].strand == "-"
        cumulative = 0
        for feature in sorted(cds_records, key=lambda item: item.start, reverse=reverse):
            expected = (3 - cumulative % 3) % 3
            if feature.phase != str(expected):
                failures.append(f"{model_id}: incorrect phase remains")
            cumulative += feature.end - feature.start + 1
    if failures:
        raise RuntimeError("; ".join(failures))

    args.output_dir.mkdir(parents=True)
    gff_output = args.output_dir / "Populus_bZIP_25_retained_candidate_models.corrected_v2.gff3"
    with gff_output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("##gff-version 3\n")
        for feature in corrected:
            handle.write(feature.to_line() + "\n")

    cds_output = args.output_dir / args.cds_fasta.name
    protein_output = args.output_dir / args.protein_fasta.name
    audit_output = args.output_dir / args.audit_tsv.name
    shutil.copy2(args.cds_fasta, cds_output)
    shutil.copy2(args.protein_fasta, protein_output)
    shutil.copy2(args.audit_tsv, audit_output)
    shutil.copy2(args.seqid_map, args.output_dir / "seqid_alias_mapping.tsv")

    write_tsv(
        args.output_dir / "gff3_phase_corrections.tsv",
        phase_changes,
        ["model_id", "mRNA_id", "CDS_index_transcript_order", "start", "end", "strand", "old_phase", "new_phase"],
    )
    write_tsv(
        args.output_dir / "gff3_added_exons.tsv",
        added_exon_rows,
        ["model_id", "mRNA_id", "exon_id", "seqid", "start", "end", "strand"],
    )
    write_tsv(
        args.output_dir / "gff3_seqid_corrections.tsv",
        seqid_changes,
        ["model_id", "feature_id", "feature_type", "old_seqid", "new_seqid"],
    )

    per_species = args.output_dir / "supported_models_by_species"
    per_species.mkdir()
    ids_by_species: dict[str, list[str]] = defaultdict(list)
    for model_id in sorted(retained_ids):
        ids_by_species[species_by_model[model_id]].append(model_id)
    for species, identifiers in sorted(ids_by_species.items()):
        species_dir = per_species / species
        species_dir.mkdir()
        with (species_dir / f"{species}.retained_candidate_models.corrected_v2.gff3").open(
            "w", encoding="utf-8", newline="\n"
        ) as handle:
            handle.write("##gff-version 3\n")
            for feature in corrected:
                if feature.model_id in identifiers:
                    handle.write(feature.to_line() + "\n")
        write_fasta(
            species_dir / f"{species}.retained_candidate_models.cds.fasta",
            identifiers,
            cds_sequences,
        )
        write_fasta(
            species_dir / f"{species}.retained_candidate_models.protein.fasta",
            identifiers,
            protein_sequences,
        )

    report = {
        "source_files_preserved": True,
        "retained_models": len(retained_ids),
        "independent_new_loci": sum(row["independent_new_gene"] == "Yes" for row in retained_rows),
        "annotation_corrections": sum(row["independent_new_gene"] == "No" for row in retained_rows),
        "gff3_counts": {
            key: sum(feature.feature_type == key for feature in corrected)
            for key in ["gene", "mRNA", "exon", "CDS"]
        },
        "phase_fields_changed": len(phase_changes),
        "phase_models_changed": len({row["model_id"] for row in phase_changes}),
        "exon_features_added": len(added_exon_rows),
        "models_receiving_exons": len({row["model_id"] for row in added_exon_rows}),
        "seqid_feature_fields_changed": len(seqid_changes),
        "seqid_models_changed": len({row["model_id"] for row in seqid_changes}),
        "cds_fasta_unchanged": sha256(cds_output) == sha256(args.cds_fasta),
        "protein_fasta_unchanged": sha256(protein_output) == sha256(args.protein_fasta),
        "internal_validation_passed": True,
    }
    (args.output_dir / "GFF3_CORRECTION_REPORT.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    files = sorted(
        path
        for path in args.output_dir.rglob("*")
        if path.is_file() and path.name != "SHA256SUMS.tsv"
    )
    with (args.output_dir / "SHA256SUMS.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("sha256\trelative_path\n")
        for path in files:
            handle.write(f"{sha256(path)}\t{path.relative_to(args.output_dir)}\n")

    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
