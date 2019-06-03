log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

condition <- snakemake@params[["contrast"]][1]
levelA <- snakemake@params[["contrast"]][2]
levelB <- snakemake@params[["contrast"]][3]

cat("Condition:",condition,"\n")
cat("Level A:",levelA,"\n")
cat("Level B:",levelB,"\n")

cat("ResultNames:",resultsNames(dds),"\n")

contrast <- c(condition,levelA,levelB)

res <- results(dds, contrast=contrast, parallel=parallel)

# shrink fold changes for lowly expressed genes #TODO Enable Shrinkage for Interactions?
res <- lfcShrink(dds, contrast=contrast, res=res)

# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=snakemake@output[["table"]])
