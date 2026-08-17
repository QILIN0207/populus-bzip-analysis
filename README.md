# Populus bZIP analysis code (v1.1.0)

This package contains the final custom scripts used for the eight analysis stages below. Figure drawing, Office formatting, temporary audits, logs and superseded scripts are not included.

| Directory | Scope |
|---|---|
| `01_candidate_models` | Candidate-model structural audit, recovery and GFF3/CDS/protein validation |
| `02_orthogroups` | OrthoFinder inputs, de novo OGG inference and sensitivity analysis |
| `03_pangene` | Pangene naming and OGG representative selection |
| `04_phylogeny_motifs` | Phylogeny, MEME/MAST and group-D heptad periodicity |
| `05_duplication_synteny` | MCScanX duplication-mode and synteny statistics |
| `06_transposable_elements` | TE overlap, matched background and effect sizes |
| `07_kaks` | Ka/Ks inputs, collision-safe IDs, YN00 execution and filtering |
| `08_rnaseq_enrichment` | RNA-seq counting, DESeq2, summaries and GO/KEGG enrichment |

## Installation

Python 3.10 or later is required for the Python scripts:

```bash
python -m pip install -r requirements.txt
```

External programs are installed separately: BITACORA 1.4, miniprot 0.18-r281, HMMER 3.4, OrthoFinder 2.5.5, DIAMOND 2.2.2, MAFFT 7.505, IQ-TREE 2.0.7, MEME Suite 5.5.9, BLAST+, MCScanX, ParaAT 2.0, KaKs_Calculator 1.2, fastp 0.23.4, HISAT2 2.2.2, SAMtools 1.23.1, HTSeq 2.1.2 and DESeq2 1.50.2.

The retained Bioconda builds:

- `hisat2-2.2.2-h503566f_0`
- `samtools-1.23.1-ha83d96e_0`
- `diamond-2.2.2-he361c42_0`

## Running the analyses

Every Python script provides `--help`. Shell wrappers print a usage line when called without arguments and refuse to overwrite an existing output directory. The principal entry points are:

- `02_orthogroups/run_orthofinder_sets.sh`
- `04_phylogeny_motifs/run_phylogeny.sh`
- `04_phylogeny_motifs/run_meme_mast.sh`
- `05_duplication_synteny/run_mcscanx.sh`
- `06_transposable_elements/te_interval_overlap.py` and `te_background_analysis.py`
- `07_kaks/run_paraat_kaks.sh`
- `08_rnaseq_enrichment/run_hisat2_htseq.sh`, `merge_htseq_counts.py`, `run_deseq2.R` and `hypergeometric_enrichment.py`

Exact inputs and outputs are listed in `SCRIPT_MANIFEST.tsv`. The small author-generated input/result files retained from v1.0.0 remain under `data/`. Third-party genomes, annotations, RNA-seq reads and full external-tool outputs are not redistributed.

## Key fixed settings

- OrthoFinder: three input sets (audited full set, no independent new loci, and annotation only).
- Phylogeny: MAFFT `--auto`; columns with >80% gaps removed; IQ-TREE ModelFinder, 1,000 UFBoot and 1,000 SH-aLRT replicates.
- MEME/MAST: protein, ZOOPS, 10 motifs, width 6-50 aa, sequence E-value 10 and motif-hit P-value 1e-4.
- MCScanX: BLASTP E-value 1e-5, five non-self hits per query, match score 50, gap penalty -1 and MATCH_SIZE 5.
- Duplication association: 1,000,000 fixed-margin Monte Carlo tables, seed 20260721, with Cramer's V.
- Ka/Ks: ParaAT/MAFFT codon alignments and YN00 in KaKs_Calculator; undefined values and the prespecified Ks/Ka/KaKs filters are recorded.
- RNA-seq: HISAT2 `--dta`; MAPQ >=15; HTSeq union mode, unstranded, exon/gene_id; DESeq2 significance requires |log2FC| >=1 and adjusted P <=0.05.

Scripts are under the MIT License; author-generated data are under CC BY 4.0.
