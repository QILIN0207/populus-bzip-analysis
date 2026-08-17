#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: run_deseq2.R COUNTS.tsv SAMPLE_METADATA.tsv CONTRASTS.tsv OUTPUT_DIR")
}
counts_path <- args[[1]]
metadata_path <- args[[2]]
contrasts_path <- args[[3]]
output_dir <- args[[4]]
if (dir.exists(output_dir)) stop("Refusing to overwrite output directory: ", output_dir)
dir.create(output_dir, recursive = TRUE)

counts <- read.delim(counts_path, check.names = FALSE, row.names = 1)
metadata <- read.delim(metadata_path, check.names = FALSE, row.names = 1)
contrasts <- read.delim(contrasts_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(colnames(counts) %in% rownames(metadata))) stop("Counts and metadata samples differ")
metadata <- metadata[colnames(counts), , drop = FALSE]
metadata$condition <- factor(metadata$condition)

counts <- counts[rowSums(counts) >= 10, , drop = FALSE]
dds <- DESeqDataSetFromMatrix(round(as.matrix(counts)), metadata, design = ~ condition)
dds <- DESeq(dds)
for (i in seq_len(nrow(contrasts))) {
  label <- contrasts$contrast_id[[i]]
  numerator <- contrasts$numerator[[i]]
  denominator <- contrasts$denominator[[i]]
  result <- results(dds, contrast = c("condition", numerator, denominator), alpha = 0.05)
  table <- as.data.frame(result)
  table$gene_id <- rownames(table)
  table$significant <- !is.na(table$padj) & table$padj <= 0.05 & abs(table$log2FoldChange) >= 1
  write.table(table, file.path(output_dir, paste0(label, ".tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
}



