log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])


svg(snakemake@output[[1]])
plotDispEsts(dds)
dev.off()

