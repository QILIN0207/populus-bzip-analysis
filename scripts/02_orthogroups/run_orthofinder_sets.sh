#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 INPUT_ROOT OUTPUT_ROOT [THREADS] [ALGORITHM_THREADS] [ORTHOFINDER]" >&2
  exit 2
}

[[ $# -ge 2 ]] || usage
input_root="$1"
output_root="$2"
threads="${3:-12}"
algorithm_threads="${4:-1}"
orthofinder="${5:-orthofinder}"

datasets=(revised_audited_25 sensitivity_no_independent_new sensitivity_annotated_only)
[[ ! -e "$output_root" ]] || { echo "Refusing to overwrite: $output_root" >&2; exit 3; }
mkdir -p "$output_root/logs" "$output_root/runs"

printf 'dataset\tinput_directory\tresult_directory\n' > "$output_root/run_manifest.tsv"
for dataset in "${datasets[@]}"; do
  input="$input_root/$dataset"
  result="$output_root/runs/$dataset"
  [[ -d "$input" ]] || { echo "Missing input directory: $input" >&2; exit 4; }
  printf '%s\t%s\t%s\n' "$dataset" "$input" "$result" >> "$output_root/run_manifest.tsv"
  "$orthofinder" -f "$input" -t "$threads" -a "$algorithm_threads" -o "$result" \
    2>&1 | tee "$output_root/logs/${dataset}.log"
  grep -q 'Done orthogroups' "$output_root/logs/${dataset}.log"
done


