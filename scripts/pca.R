library("DESeq2")

# load deseq2 data
dds <- load(snakemake@input[[1]])

# obtain normalized counts
ntd <- normTransform(dds)
pdf(snakemake@output[[1]])
plotPCA(ntd, intgroup=snakemake@params[["pca_labels"]])
dev.off()
