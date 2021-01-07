#! /usr/bin/env Rscript
library(DESeq2)

args <- commandArgs()
gene_reads_count <- args[6]
trans_reads_count <- args[7]
mata_file <- args[8]
gene_diff_res <- args[9]
trans_diff_res <- args[10]
#print(gene_reads_count)
#print(trans_reads_count)
#print(mata_file)
#print(gene_diff_res)
#print(trans_diff_res)

meta_dt <- read.csv(mata_file, head=TRUE, row.names="sample")

gene_dt <- read.csv(gene_reads_count, head=TRUE, row.names="gene_id")
gene_dds <- DESeqDataSetFromMatrix(countData=gene_dt, colData=meta_dt, design=~group)
gene_dds <- DESeq(gene_dds)
gene_res <- results(gene_dds)
gene_res <- gene_res[order(gene_res$padj), ]
write.csv(gene_res, file=gene_diff_res)

trans_dt <- read.csv(trans_reads_count, head=TRUE, row.names="transcript_id")
trans_dds <- DESeqDataSetFromMatrix(countData=trans_dt, colData=meta_dt, design=~group)
trans_dds <- DESeq(trans_dds)
trans_res <- results(trans_dds)
trans_res <- trans_res[order(trans_res$padj), ]
write.csv(trans_res, file=trans_diff_res)
