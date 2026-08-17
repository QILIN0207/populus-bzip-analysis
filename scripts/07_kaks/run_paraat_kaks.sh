#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 PROTEIN_FASTA CDS_FASTA PAIRS_TSV OUTPUT_DIR PARAAT_PL KAKS_CALCULATOR [THREADS] [MAFFT]" >&2
  exit 2
}

[[ $# -ge 6 ]] || usage
protein="$1"
cds="$2"
pairs="$3"
output_dir="$4"
paraat="$5"
kaks="$6"
threads="${7:-8}"
mafft="${8:-mafft}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 3; }
mkdir -p "$output_dir/logs"
python3 "$script_dir/make_paraat_safe_inputs.py" --protein "$protein" --cds "$cds" \
  --pairs "$pairs" --output-dir "$output_dir/safe_input"
printf '%s\n' "$threads" > "$output_dir/proc.txt"
perl "$paraat" -h "$output_dir/safe_input/paraat_safe.pairs.tsv" \
  -a "$output_dir/safe_input/paraat_safe.protein.fasta" \
  -n "$output_dir/safe_input/paraat_safe.cds.fasta" -p "$output_dir/proc.txt" \
  -o "$output_dir/paraat_axt" -m "$mafft" -f axt > "$output_dir/logs/paraat.log" 2>&1
grep -q 'Mission Accomplished' "$output_dir/logs/paraat.log"

mkdir -p "$output_dir/kaks_alias_results"
export kaks output_dir
find "$output_dir/paraat_axt" -maxdepth 1 -type f -name '*.axt' -print0 | sort -z | \
  xargs -0 -P "$threads" -I '{}' bash -c '
    input="$1"; base="$(basename "$input" .axt)"
    "$kaks" -i "$input" -o "$output_dir/kaks_alias_results/${base}.kaks.tsv" -m YN
  ' _ '{}'

