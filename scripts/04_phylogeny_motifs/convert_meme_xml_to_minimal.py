#!/usr/bin/env python3
"""Convert MEME XML output to validated Minimal MEME Motif Format."""

from __future__ import annotations

import argparse
import math
import xml.etree.ElementTree as ET
from pathlib import Path


PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-xml", type=Path, required=True)
    parser.add_argument("--output-meme", type=Path, required=True)
    return parser.parse_args()


def values_by_letter(element: ET.Element) -> dict[str, str]:
    values: dict[str, str] = {}
    for value in element.findall("value"):
        letter = value.attrib["letter_id"]
        if value.text is None:
            raise ValueError(f"Missing value for letter {letter}")
        values[letter] = value.text.strip()
    return values


def validate_probability_row(
    values: list[str], label: str, *, tolerance: float = 5e-5
) -> None:
    numbers = [float(value) for value in values]
    if any(not math.isfinite(value) or value < 0 or value > 1 for value in numbers):
        raise ValueError(f"Invalid probability in {label}")
    if not math.isclose(sum(numbers), 1.0, rel_tol=0, abs_tol=tolerance):
        raise ValueError(f"Probabilities do not sum to 1 in {label}: {sum(numbers)}")


def main() -> None:
    args = parse_args()
    if not args.input_xml.is_file():
        raise FileNotFoundError(args.input_xml)
    if args.output_meme.exists():
        raise FileExistsError(f"Refusing overwrite: {args.output_meme}")

    root = ET.parse(args.input_xml).getroot()
    if root.tag != "MEME":
        raise ValueError(f"Expected MEME XML root, found {root.tag}")
    version = root.attrib.get("version", "5")

    training_set = root.find("training_set")
    if training_set is None:
        raise ValueError("Missing training_set")
    alphabet = training_set.find("alphabet")
    if alphabet is None or alphabet.attrib.get("name", "").lower() != "protein":
        raise ValueError("Expected a protein MEME result")
    defined_letters = {letter.attrib["id"] for letter in alphabet.findall("letter")}
    missing_letters = set(PROTEIN_ALPHABET) - defined_letters
    if missing_letters:
        raise ValueError(f"Missing protein letters: {sorted(missing_letters)}")

    background_element = root.find("./model/background_frequencies/alphabet_array")
    if background_element is None:
        raise ValueError("Missing background frequencies")
    background = values_by_letter(background_element)
    background_values = [background[letter] for letter in PROTEIN_ALPHABET]
    validate_probability_row(background_values, "background", tolerance=1e-3)
    background_total = sum(float(value) for value in background_values)
    normalized_background = {
        letter: float(background[letter]) / background_total
        for letter in PROTEIN_ALPHABET
    }

    motif_elements = root.findall("./motifs/motif")
    if not motif_elements:
        raise ValueError("No motifs found")

    lines = [
        f"MEME version {version}",
        "",
        f"ALPHABET= {PROTEIN_ALPHABET}",
        "",
        f"Background letter frequencies (from {args.input_xml.name}):",
    ]
    lines.append(
        " ".join(
            f"{letter} {normalized_background[letter]:.10f}"
            for letter in PROTEIN_ALPHABET
        )
    )
    lines.append("")

    for index, motif in enumerate(motif_elements, start=1):
        motif_id = motif.attrib.get("id", f"motif_{index}")
        alternate = motif.attrib.get("alt", f"MEME-{index}")
        width = int(motif.attrib["width"])
        sites = motif.attrib.get("sites", "20")
        e_value = motif.attrib.get("e_value", "0")
        matrix = motif.findall("./probabilities/alphabet_matrix/alphabet_array")
        if len(matrix) != width:
            raise ValueError(
                f"Motif {motif_id} declares width {width} but has {len(matrix)} rows"
            )
        lines.extend(
            [
                f"MOTIF {motif_id} {alternate}",
                (
                    "letter-probability matrix: "
                    f"alength= 20 w= {width} nsites= {sites} E= {e_value}"
                ),
            ]
        )
        for row_index, row in enumerate(matrix, start=1):
            row_values = values_by_letter(row)
            values = [row_values[letter] for letter in PROTEIN_ALPHABET]
            validate_probability_row(values, f"{motif_id} row {row_index}")
            lines.append(" ".join(values))
        lines.append("")

    args.output_meme.parent.mkdir(parents=True, exist_ok=True)
    args.output_meme.write_text("\n".join(lines) + "\n", encoding="ascii")
    print(
        f"Wrote {len(motif_elements)} protein motifs to {args.output_meme} "
        f"using MEME version {version}."
    )


if __name__ == "__main__":
    main()
