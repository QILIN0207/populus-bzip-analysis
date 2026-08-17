#!/usr/bin/env python3
"""Prepare revised Populus bZIP protein inputs for MEME and MAST uploads."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path


SPECIES_CODES = {
    "Populus_adenopoda": "Pade",
    "Populus_alba": "Palb",
    "Populus_alba_var__pyramidalis": "Pavp",
    "Populus_davidiana": "Pdav",
    "Populus_deltoide": "Pdel",
    "Populus_deltoides": "Pdel",
    "Populus_euphratica": "Peup",
    "Populus_ilicifolia": "Pili",
    "Populus_koreana": "Pkor",
    "Populus_lasiocarpa": "Plas",
    "Populus_pruinosa": "Ppru",
    "Populus_pseudoglauca": "Ppsg",
    "Populus_qiongdaoensis": "Pqio",
    "Populus_rotundifolia": "Prot",
    "Populus_simonii": "Psim",
    "Populus_szechuanica": "Psze",
    "Populus_tremula": "Ptre",
    "Populus_trichocarpa": "Ptri",
    "Populus_wuana": "Pwua",
    "Populus_yunnanensis": "Pyun",
}

ALLOWED_PROTEIN = set("ACDEFGHIKLMNPQRSTVWYX")


@dataclass(frozen=True)
class FastaRecord:
    identifier: str
    header: str
    sequence: str
    species: str = ""
    source_file: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--representatives", type=Path, required=True)
    parser.add_argument("--primary-dir", type=Path, required=True)
    parser.add_argument("--no-new-dir", type=Path, required=True)
    parser.add_argument("--annotated-only-dir", type=Path, required=True)
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--unassigned", type=Path, required=True)
    parser.add_argument("--nomenclature", type=Path, required=True)
    parser.add_argument("--subfamilies", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def read_fasta(path: Path, species: str = "") -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    chunks: list[str] = []

    def finish() -> None:
        nonlocal header, chunks
        if header is None:
            return
        sequence = "".join(chunks).replace(" ", "").upper()
        identifier = header.split()[0]
        if not sequence:
            raise ValueError(f"Empty sequence for {identifier} in {path}")
        invalid = sorted(set(sequence) - ALLOWED_PROTEIN)
        if invalid:
            raise ValueError(
                f"Unsupported protein characters for {identifier} in {path}: {invalid}"
            )
        records.append(
            FastaRecord(identifier, header, sequence, species, path.name)
        )

    with path.open(encoding="utf-8-sig") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                finish()
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before first header in {path}")
                chunks.append(line)
    finish()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    identifiers = [record.identifier for record in records]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError(f"Duplicate FASTA IDs in {path}")
    return records


def read_dataset(directory: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    observed: set[str] = set()
    for species in SPECIES_CODES:
        path = directory / f"{species}.fa"
        if not path.is_file():
            raise FileNotFoundError(path)
        species_records = read_fasta(path, species)
        expected_prefixes = {species, species.replace("_var__", "_var_")}
        for record in species_records:
            if record.identifier in observed:
                raise ValueError(f"Duplicate dataset ID: {record.identifier}")
            if not any(
                record.identifier.startswith(f"{prefix}__")
                for prefix in expected_prefixes
            ):
                raise ValueError(
                    f"ID {record.identifier} does not match species file {path.name}"
                )
            observed.add(record.identifier)
            records.append(record)
    extras = sorted(
        path.name
        for path in directory.glob("*.fa")
        if path.stem not in SPECIES_CODES
    )
    if extras:
        raise ValueError(f"Unexpected FASTA files in {directory}: {extras}")
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_gene_to_orthogroup(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    rows = read_tsv(path)
    for row in rows:
        orthogroup = row["Orthogroup"]
        for species, cell in row.items():
            if species == "Orthogroup" or not cell:
                continue
            for identifier in (item.strip() for item in cell.split(",")):
                if not identifier:
                    continue
                if identifier in mapping:
                    raise ValueError(f"Gene assigned more than once: {identifier}")
                mapping[identifier] = orthogroup
    return mapping


def write_fasta(path: Path, records: list[tuple[str, str]], width: int = 60) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for identifier, sequence in records:
            handle.write(f">{identifier}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    if args.output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists; refusing overwrite: {args.output_dir}"
        )
    args.output_dir.mkdir(parents=True)

    primary = read_dataset(args.primary_dir)
    no_new = read_dataset(args.no_new_dir)
    annotated_only = read_dataset(args.annotated_only_dir)
    expected_counts = {
        "primary": 1764,
        "no_independent_new": 1743,
        "annotated_only": 1739,
    }
    observed_counts = {
        "primary": len(primary),
        "no_independent_new": len(no_new),
        "annotated_only": len(annotated_only),
    }
    if observed_counts != expected_counts:
        raise ValueError(
            f"Unexpected dataset counts: {observed_counts}; expected {expected_counts}"
        )

    primary_by_id = {record.identifier: record for record in primary}
    no_new_by_id = {record.identifier: record for record in no_new}
    annotated_by_id = {record.identifier: record for record in annotated_only}
    for label, subset in (
        ("no_independent_new", no_new_by_id),
        ("annotated_only", annotated_by_id),
    ):
        missing = sorted(set(subset) - set(primary_by_id))
        if missing:
            raise ValueError(f"{label} contains IDs absent from primary: {missing}")
        sequence_mismatches = sorted(
            identifier
            for identifier, record in subset.items()
            if record.sequence != primary_by_id[identifier].sequence
        )
        if sequence_mismatches:
            raise ValueError(
                f"{label} sequence mismatches against primary: {sequence_mismatches}"
            )

    nomenclature_rows = read_tsv(args.nomenclature)
    og_to_nomenclature = {
        row["New_Orthogroup"]: row for row in nomenclature_rows
    }
    if len(og_to_nomenclature) != 79:
        raise ValueError("Expected 79 revised orthogroups in nomenclature table")
    subfamily_rows = read_tsv(args.subfamilies)
    ptbzip_to_subfamily = {row["PtbZIP_ID"]: row for row in subfamily_rows}
    if len(ptbzip_to_subfamily) != 79:
        raise ValueError("Expected 79 PtbZIP rows in subfamily table")

    gene_to_og = build_gene_to_orthogroup(args.orthogroups)
    unassigned_to_og = build_gene_to_orthogroup(args.unassigned)
    overlap = set(gene_to_og) & set(unassigned_to_og)
    if overlap:
        raise ValueError(f"Genes both assigned and unassigned: {sorted(overlap)}")
    if len(gene_to_og) != 1762 or len(unassigned_to_og) != 2:
        raise ValueError(
            f"Unexpected assignment counts: assigned={len(gene_to_og)}, "
            f"unassigned={len(unassigned_to_og)}"
        )
    if set(primary_by_id) != set(gene_to_og) | set(unassigned_to_og):
        missing = sorted(set(primary_by_id) - set(gene_to_og) - set(unassigned_to_og))
        extra = sorted((set(gene_to_og) | set(unassigned_to_og)) - set(primary_by_id))
        raise ValueError(f"Orthogroup mapping mismatch; missing={missing}, extra={extra}")

    upload_ids: dict[str, str] = {}
    species_counters = {species: 0 for species in SPECIES_CODES}
    for record in primary:
        species_counters[record.species] += 1
        upload_ids[record.identifier] = (
            f"{SPECIES_CODES[record.species]}R_bZIP"
            f"{species_counters[record.species]:03d}"
        )
    if len(upload_ids) != len(set(upload_ids.values())):
        raise ValueError("Generated upload IDs are not unique")

    primary_fasta = args.output_dir / "02_MAST_New1764_primary_audited_proteins.fasta"
    no_new_fasta = (
        args.output_dir / "03_MAST_New1743_no_independent_new_proteins.fasta"
    )
    annotated_fasta = (
        args.output_dir / "04_MAST_New1739_annotated_only_proteins.fasta"
    )
    write_fasta(
        primary_fasta,
        [(upload_ids[r.identifier], r.sequence) for r in primary],
    )
    write_fasta(
        no_new_fasta,
        [
            (upload_ids[r.identifier], r.sequence)
            for r in primary
            if r.identifier in no_new_by_id
        ],
    )
    write_fasta(
        annotated_fasta,
        [
            (upload_ids[r.identifier], r.sequence)
            for r in primary
            if r.identifier in annotated_by_id
        ],
    )

    gene_mapping_rows: list[dict[str, object]] = []
    for record in primary:
        identifier = record.identifier
        if identifier in gene_to_og:
            orthogroup = gene_to_og[identifier]
            nomenclature = og_to_nomenclature.get(orthogroup)
            if nomenclature is None:
                raise ValueError(f"No PtbZIP nomenclature for {orthogroup}")
            ptbzip = nomenclature["Pangene_ID"]
            subfamily = ptbzip_to_subfamily[ptbzip]
            assignment_status = "assigned"
        else:
            orthogroup = unassigned_to_og[identifier]
            nomenclature = {}
            ptbzip = ""
            subfamily = {}
            assignment_status = "unassigned"
        gene_mapping_rows.append(
            {
                "upload_id": upload_ids[identifier],
                "original_id": identifier,
                "species": record.species,
                "species_code": SPECIES_CODES[record.species],
                "source_file": record.source_file,
                "sequence_length": len(record.sequence),
                "new_orthogroup": orthogroup,
                "PtbZIP_ID": ptbzip,
                "pan_genome_class": nomenclature.get("Category", ""),
                "subfamily": subfamily.get("subfamily", ""),
                "subfamily_color": subfamily.get("subfamily_color", ""),
                "assignment_status": assignment_status,
                "included_primary_1764": "yes",
                "included_no_independent_new_1743": (
                    "yes" if identifier in no_new_by_id else "no"
                ),
                "included_annotated_only_1739": (
                    "yes" if identifier in annotated_by_id else "no"
                ),
                "X_residue_count": record.sequence.count("X"),
            }
        )
    gene_mapping = args.output_dir / "02_MAST_New1764_gene_id_mapping.tsv"
    write_tsv(
        gene_mapping,
        list(gene_mapping_rows[0]),
        gene_mapping_rows,
    )

    representative_records = read_fasta(args.representatives)
    if len(representative_records) != 79:
        raise ValueError(
            f"Expected 79 representative sequences, found {len(representative_records)}"
        )
    representative_fasta_rows: list[tuple[str, str]] = []
    representative_mapping_rows: list[dict[str, object]] = []
    seen_pangenes: set[str] = set()
    for record in representative_records:
        pangene = record.identifier.split("__", 1)[0]
        if pangene in seen_pangenes:
            raise ValueError(f"Duplicate representative pangene: {pangene}")
        seen_pangenes.add(pangene)
        if pangene not in ptbzip_to_subfamily:
            raise ValueError(f"Representative pangene missing subfamily: {pangene}")
        if "__" not in record.identifier:
            raise ValueError(f"Representative lacks member ID: {record.identifier}")
        member_id = record.identifier.split("__", 1)[1]
        if member_id not in primary_by_id:
            raise ValueError(f"Representative member absent from primary input: {member_id}")
        if primary_by_id[member_id].sequence != record.sequence:
            raise ValueError(f"Representative sequence mismatch: {pangene}")
        representative_fasta_rows.append((pangene, record.sequence))
        subfamily = ptbzip_to_subfamily[pangene]
        representative_mapping_rows.append(
            {
                "PtbZIP_ID": pangene,
                "original_representative_header": record.identifier,
                "representative_member_id": member_id,
                "representative_species": primary_by_id[member_id].species,
                "sequence_length": len(record.sequence),
                "subfamily": subfamily["subfamily"],
                "subfamily_color": subfamily["subfamily_color"],
                "pan_genome_class": subfamily["pan_genome_class"],
                "X_residue_count": record.sequence.count("X"),
            }
        )
    if seen_pangenes != set(ptbzip_to_subfamily):
        raise ValueError("Representative and subfamily PtbZIP sets differ")

    representative_fasta = (
        args.output_dir / "01_MEME_New79_PtbZIP_representatives_protein.fasta"
    )
    representative_mapping = (
        args.output_dir / "01_MEME_New79_representatives_mapping.tsv"
    )
    write_fasta(representative_fasta, representative_fasta_rows)
    write_tsv(
        representative_mapping,
        list(representative_mapping_rows[0]),
        representative_mapping_rows,
    )

    readme = args.output_dir / "README_UPLOAD_INSTRUCTIONS.txt"
    readme.write_text(
        """Revised Populus bZIP MEME/MAST upload inputs

UPLOAD NOW TO MEME
  File: 01_MEME_New79_PtbZIP_representatives_protein.fasta
  Alphabet: Protein
  Site distribution: ZOOPS (zero or one occurrence per sequence)
  Maximum number of motifs: 10
  Minimum motif width: 6 aa
  Maximum motif width: 50 aa
  Objective function: Classic
  Background: 0-order frequencies calculated from the input sequences
  E-value stopping threshold: no explicit threshold; request exactly 10 motifs

DOWNLOAD AND RETAIN FROM MEME
  Download the complete result archive and retain meme.xml, meme.txt,
  meme.html, logos, submitted parameters, MEME version, and job identifier.

USE LATER WITH THE NEW meme.xml
  02_MAST_New1764_primary_audited_proteins.fasta
  03_MAST_New1743_no_independent_new_proteins.fasta
  04_MAST_New1739_annotated_only_proteins.fasta

The three MAST files share stable upload IDs. Use
02_MAST_New1764_gene_id_mapping.tsv to restore original gene IDs, revised
PtbZIP IDs, pan-genome classes, and subfamilies. The primary file contains
1,764 proteins. Subfamily-frequency summaries use the 1,762 assigned proteins;
the two unassigned proteins are retained in the scan and reported separately.

Do not upload the 1,764-protein file for de novo motif discovery. It is the
target database for scanning the motifs discovered from the 79 representatives.
""",
        encoding="utf-8",
        newline="\n",
    )
    output_files = sorted(
        path
        for path in args.output_dir.iterdir()
        if path.is_file() and path.name != "SHA256SUMS.tsv"
    )
    qa = {
        "all_passed": True,
        "representatives": {
            "sequences": len(representative_fasta_rows),
            "unique_ids": len(seen_pangenes),
            "total_residues": sum(len(seq) for _, seq in representative_fasta_rows),
            "X_residues": sum(seq.count("X") for _, seq in representative_fasta_rows),
        },
        "primary": {
            "sequences": len(primary),
            "assigned": len(gene_to_og),
            "unassigned": len(unassigned_to_og),
            "total_residues": sum(len(record.sequence) for record in primary),
            "X_residues": sum(record.sequence.count("X") for record in primary),
        },
        "sensitivity_no_independent_new": {"sequences": len(no_new)},
        "sensitivity_annotated_only": {"sequences": len(annotated_only)},
        "species_counts": species_counters,
        "unassigned_original_ids": sorted(unassigned_to_og),
        "protein_alphabet": "ACDEFGHIKLMNPQRSTVWYX",
    }
    qa_path = args.output_dir / "QA_REPORT.json"
    qa_path.write_text(
        json.dumps(qa, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    output_files.append(qa_path)

    checksums = args.output_dir / "SHA256SUMS.tsv"
    with checksums.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("sha256\tfile\n")
        for path in sorted(output_files):
            handle.write(f"{sha256(path)}\t{path.name}\n")

    print(json.dumps(qa, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
