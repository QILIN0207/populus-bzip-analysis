#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 PREFIX PROTEIN_FASTA OUTPUT_DIR [THREADS] [BLASTP] [MCScanX] [CLASSIFIER]" >&2
  echo "PREFIX.gff must be a four-column MCScanX coordinate file." >&2
  exit 2
}

[[ $# -ge 3 ]] || usage
prefix="$1"
protein="$2"
output_dir="$3"
threads="${4:-8}"
blastp="${5:-blastp}"
mcscanx="${6:-MCScanX}"
classifier="${7:-duplicate_gene_classifier}"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

[[ -s "${prefix}.gff" && -s "$protein" ]] || { echo "Missing ${prefix}.gff or protein FASTA" >&2; exit 3; }
[[ ! -e "$output_dir" ]] || { echo "Refusing to overwrite: $output_dir" >&2; exit 4; }
mkdir -p "$output_dir"
cp "${prefix}.gff" "$output_dir/analysis.gff"

makeblastdb -in "$protein" -dbtype prot -out "$output_dir/proteins" >/dev/null
"$blastp" -query "$protein" -db "$output_dir/proteins" -outfmt 6 -evalue 1e-5 \
  -num_threads "$threads" -max_target_seqs 100000 -out "$output_dir/analysis.blast.raw"
python3 "$script_dir/filter_mcscanx_blast.py" "$output_dir/analysis.blast.raw" \
  "$output_dir/analysis.blast" --max-hits 5 --evalue 1e-5
(cd "$output_dir" && "$mcscanx" analysis -k 50 -g -1 -s 5 -m 25 -e 1e-5 -b 2)
(cd "$output_dir" && "$classifier" analysis)

