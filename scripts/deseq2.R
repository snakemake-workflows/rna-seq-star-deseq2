library("DESeq2")

dds <- load(snakemake@input[[1]])

contrast <- c("condition", snakemake@params[["contrast"]])
res <- results(dds, contrast=contrast)
# shrink fold changes for lowly expressed genes
res <- lfcShrink(dds, contrast=contrast, res=res)
# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
pdf(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=snakemake@output[["table"]])
