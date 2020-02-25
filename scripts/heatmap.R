log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("pheatmap")
library("vsn")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

ntd <- normTransform(dds)

coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c(colnames(coldata))])
svg(snakemake@output[[1]])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)
dev.off()

