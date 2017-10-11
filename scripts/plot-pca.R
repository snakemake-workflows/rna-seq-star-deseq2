sink(file=snakemake@log[[1]], type="message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

# obtain normalized counts
counts <- rlog(dds, blind=FALSE)
svg(snakemake@output[[1]])
plotPCA(counts, intgroup=snakemake@params[["pca_labels"]])
dev.off()
