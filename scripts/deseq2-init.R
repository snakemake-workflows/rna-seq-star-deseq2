#sink(file=snakemake@log[[1]])

library("DESeq2")
library("BiocParallel")

# setup parallelization
register(MulticoreParam(snakemake@threads))

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene")
coldata <- read.table(snakemake@input[["samples"]], header=TRUE, row.names="sample")

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# TODO optionally allow to collapse technical replicates
dds <- DESeq(dds, parallel=TRUE)

saveRDS(dds, file=snakemake@output[[1]])
