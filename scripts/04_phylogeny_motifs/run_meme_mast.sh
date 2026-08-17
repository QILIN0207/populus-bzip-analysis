#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 REPRESENTATIVES_FASTA ALL_PROTEINS_FASTA OUTPUT_DIR [THREADS] [MEME] [MAST]" >&2
  exit 2
}

[[ $# -ge 3 ]] || usage
representatives="$1"
all_proteins="$2"
output_dir="$3"
threads="${4:-4}"
meme_bin="${5:-meme}"
mast_bin="${6:-mast}"

[[ -s "$representatives" && -s "$all_proteins" ]] || { echo "Missing FASTA input" >&2; exit 3; }
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 4; }
mkdir -p "$output_dir"

"$meme_bin" "$representatives" -protein -oc "$output_dir/meme" -mod zoops \
  -nmotifs 10 -minw 6 -maxw 50 -p "$threads"
"$mast_bin" "$output_dir/meme/meme.txt" "$all_proteins" -oc "$output_dir/mast" \
  -ev 10 -mt 0.0001


