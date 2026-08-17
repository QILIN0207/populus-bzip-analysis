#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 INPUT_FASTA OUTPUT_DIR [THREADS] [MAFFT] [IQTREE]" >&2
  exit 2
}

[[ $# -ge 2 ]] || usage
input_fasta="$1"
output_dir="$2"
threads="${3:-6}"
mafft="${4:-mafft}"
iqtree="${5:-iqtree2}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ -s "$input_fasta" ]] || { echo "Missing input FASTA: $input_fasta" >&2; exit 3; }
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 4; }
mkdir -p "$output_dir"

alignment="$output_dir/full_length.mafft.fasta"
trimmed="$output_dir/full_length.mafft.gap80.fasta"
report="$output_dir/gap80_trim_report.tsv"

"$mafft" --auto --thread "$threads" "$input_fasta" > "$alignment"
python3 "$script_dir/trim_alignment_by_gap_fraction.py" \
  --alignment "$alignment" --output "$trimmed" --report "$report" --max-gap-fraction 0.80
"$iqtree" -s "$trimmed" -m MFP -B 1000 -alrt 1000 -T "$threads" \
  --prefix "$output_dir/bzip_full_length_gap80"

