# Populus bZIP analysis scripts and data

This repository provides the custom scripts and associated data requested for reproducing the following analyses:

1. overlap between transposable elements (TEs) and bZIP gene analysis windows, including the matched non-bZIP background comparison and OGG-level effect-size analysis;
2. one-sided hypergeometric GO and KEGG enrichment tests with Benjamini-Hochberg correction; and
3. orthogroup occupancy and copy-number statistics.

The repository also contains the 25 retained gene models and the OrthoFinder assignment files used in the analysis.

## Repository structure

```text
scripts/
  te_interval_overlap.py
  te_background_analysis.py
  te_io.py
  hypergeometric_enrichment.py
  summarize_orthogroups.py
data/
  gene_models/
  orthogroups/
  te_overlap/
  te_background/
  enrichment/
```

## Requirements

Python 3.10 or later is required. Install the Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

## TE interval overlap

```bash
python scripts/te_interval_overlap.py \
  --genes data/te_overlap/gene_windows.tsv \
  --tes data/te_overlap/te_intervals.tsv \
  --output-events results/te_overlap_events.tsv \
  --output-summary results/gene_te_summary.tsv
```

Coordinates in the two input files are zero-based and half-open. TE coverage is calculated from the union of overlapping intervals within each gene analysis window.

## Matched TE background and OGG comparison

This analysis requires the species GFF3 files, FASTA index files, sequence-identifier mappings, HiTE libraries and RepeatMasker outputs described in the manuscript. These third-party genome-scale files are not redistributed.

```bash
python scripts/te_background_analysis.py \
  --bzip-summary data/te_background/bzip_gene_summary.tsv \
  --ogg-metrics data/te_background/ogg_te_metrics.tsv \
  --te-root /path/to/species_te_directories \
  --mcscanx-root /path/to/species_gff_directories \
  --output-root results/te_background \
  --matching-replicates 1000 \
  --bootstrap-replicates 50000 \
  --seed 20260810
```

The final matched-background, per-genome TE annotation and OGG-level effect-size tables are provided in `data/te_background/`.

## GO and KEGG enrichment

```bash
python scripts/hypergeometric_enrichment.py \
  --input data/enrichment/enrichment_counts.tsv \
  --output results/enrichment_results.tsv
```

Benjamini-Hochberg correction is applied separately to GO and KEGG tests.

## Orthogroup statistics

```bash
python scripts/summarize_orthogroups.py \
  --orthogroups data/orthogroups/Orthogroups.tsv \
  --unassigned data/orthogroups/Orthogroups_UnassignedGenes.tsv \
  --output-classification results/ogg_classification.tsv \
  --output-summary results/ogg_summary.json
```

Across the 19 genomes, core OGGs occur in 19 genomes, soft-core OGGs in 17 or 18 genomes, shell OGGs in 2-16 genomes, and cloud OGGs in one genome. An OGG is classified as copy-number-variable when its member count differs among the 19 genomes.

## Data files

Detailed file definitions are provided in `data/README.md`. Original gene and transcript identifiers are preserved in the OrthoFinder files.

Third-party genome assemblies, annotations, RNA-seq reads and complete RepeatMasker outputs remain subject to their original access and licensing terms.

## Licences

Scripts are released under the MIT License. Author-generated data files are released under the Creative Commons Attribution 4.0 International License.
