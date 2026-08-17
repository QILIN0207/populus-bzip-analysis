# Data files

## gene_models

- `retained_models.gff3`: genomic features for the 25 retained models.
- `retained_models.cds.fasta`: coding sequences for the 25 retained models.
- `retained_models.protein.fasta`: protein sequences for the 25 retained models.
- `model_manifest.tsv`: standardized transcript identifier, species and model type. The set comprises 21 independent new loci and four annotation corrections.

## orthogroups

- `Orthogroups.tsv`: OrthoFinder OGG membership.
- `Orthogroups.GeneCount.tsv`: per-genome gene counts for each OGG.
- `Orthogroups_UnassignedGenes.tsv`: proteins not assigned to an OGG.
- `Orthogroup_to_PtbZIP.tsv`: OGG-to-pangene mapping used in the manuscript.

Original OrthoFinder member identifiers are retained unchanged. Suffixes present in those identifiers are therefore part of the deposited analysis input rather than repository workflow labels.

## te_overlap

- `gene_windows.tsv`: one analysis window per bZIP gene.
- `te_intervals.tsv`: unique TE intervals intersecting at least one analysis window.

Coordinates are zero-based and half-open.

## te_background

- `bzip_gene_summary.tsv`: bZIP loci used as matching targets.
- `matched_background_summary.tsv`: observed bZIP and matched non-bZIP TE-overlap statistics.
- `per_genome_te_annotation_summary.tsv`: per-genome TE annotation totals.
- `ogg_te_metrics.tsv`: OGG-level TE metrics.
- `cliffs_delta_bootstrap_ci.tsv`: Mann-Whitney U results, Cliff's delta and bootstrap confidence intervals.

## enrichment

`enrichment_counts.tsv` contains the foreground and background counts used for the one-sided GO and KEGG hypergeometric tests.
