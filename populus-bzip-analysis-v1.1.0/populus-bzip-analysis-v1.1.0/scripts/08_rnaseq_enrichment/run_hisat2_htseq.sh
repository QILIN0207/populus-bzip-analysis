#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 SAMPLE_MANIFEST.tsv OUTPUT_DIR [THREADS]" >&2
  echo "Columns: sample_id, index_prefix, gtf, read1, read2 (read2 may be NA)." >&2
  exit 2
}

[[ $# -ge 2 ]] || usage
manifest="$1"
output_dir="$2"
threads="${3:-10}"
[[ -s "$manifest" ]] || { echo "Missing manifest: $manifest" >&2; exit 3; }
mkdir -p "$output_dir"

tail -n +2 "$manifest" | while IFS=$'\t' read -r sample index_prefix gtf read1 read2; do
  [[ -n "$sample" ]] || continue
  sample_dir="$output_dir/$sample"
  [[ ! -e "$sample_dir" ]] || { echo "Refusing to overwrite sample: $sample_dir" >&2; exit 4; }
  mkdir -p "$sample_dir"
  if [[ -n "${read2:-}" && "$read2" != "NA" ]]; then
    hisat2 -p "$threads" --dta -x "$index_prefix" -1 "$read1" -2 "$read2" \
      -S "$sample_dir/alignment.sam" 2> "$sample_dir/hisat2.log"
  else
    hisat2 -p "$threads" --dta -x "$index_prefix" -U "$read1" \
      -S "$sample_dir/alignment.sam" 2> "$sample_dir/hisat2.log"
  fi
  samtools view -bS -q 15 -F 4 "$sample_dir/alignment.sam" | \
    samtools sort -o "$sample_dir/alignment.coordinate_sorted.bam"
  samtools index "$sample_dir/alignment.coordinate_sorted.bam"
  samtools sort -n -o "$sample_dir/alignment.name_sorted.bam" \
    "$sample_dir/alignment.coordinate_sorted.bam"
  htseq-count -m union -r name -s no -f bam -t exon -i gene_id \
    "$sample_dir/alignment.name_sorted.bam" "$gtf" > "$sample_dir/counts.tsv"
  rm "$sample_dir/alignment.sam"
done


