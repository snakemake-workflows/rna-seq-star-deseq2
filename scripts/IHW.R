log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("pracma")
library("DESeq2")
library("IHW")
library("ggplot2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])

IHWalpha <- as.numeric(snakemake@params[["IHWalpha"]])

res <- results(dds, parallel=parallel)
# sort by p-value
res <- res[order(res$padj),]
res <- as.data.frame(res)

# histogram 
svg(snakemake@output[["pvalHisto1"]])
ggplot(res, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
dev.off()
# IHW
ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = IHWalpha)
ggplot(as.data.frame(ihwRes), aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
dev.off()
write.table(as.data.frame(ihwRes), file=snakemake@output[["IHW_data"]])
# estimatedWeights
svg(snakemake@output[["IHW_plots"]])
plot(ihwRes)
dev.off()

# decisionBoundary
svg(snakemake@output[["IHW_plots2"]])
plot(ihwRes, what = "decisionboundary")
dev.off()

