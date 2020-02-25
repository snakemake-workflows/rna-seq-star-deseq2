log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("pheatmap")
library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])
thr <- as.numeric(snakemake@params[["thr"]])

resTC <- results(dds)
betas <- coef(dds)
topGenes <- head(order(resTC$padj),20)
# should remove the first and second column from the betas which ideally would be intercept and another comparison (might not work as intended)
mat <- betas[topGenes, -c(1,2)]
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr

svg(snakemake@output[[1]])
pheatmap(mat, breaks=seq(from=-thr, to=thr),
         cluster_col=FALSE)

while (!is.null(dev.list()))  dev.off()