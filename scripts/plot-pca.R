#sink(file=snakemake@log[[1]])

library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- counts(dds, normalized=TRUE)
svg(snakemake@output[[1]])
plotPCA(counts, intgroup=snakemake@params[["pca_labels"]])
dev.off()
