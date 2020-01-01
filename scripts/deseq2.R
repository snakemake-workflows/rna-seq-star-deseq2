log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("tidyverse")
library("apeglm")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])


# # creates a plot for each coef
# createPlots <- function(list) {    
#   for (i in 2:length(list)) { 
#     res <- lfcShrink(tryDESeq, coef=i, type="apeglm")
#     plotMA(res, main="apeglm")
#   }
# }




# creates a plot for given coef
coef <- snakemake@params[["contrast"]]
res <- lfcShrink(tryDESeq, coef=coef, type="apeglm")
plotMA(res, main="apeglm")
res <- results(dds, name=coef, parallel=parallel)
res <- res[order(res$padj),]
svg(snakemake@output[["ma_plot"]])
dev.off()
write.table(as.data.frame(res), file=snakemake@output[["table"]])


#old
#contrast <- c("condition", snakemake@params[["contrast"]])
#new, still needs to be fixxed
#contrast <- c("condition", "treated", "untreated")

# res <- results(dds, parallel=parallel)
# #res <- results(dds, contrast=contrast, parallel=parallel)
# # shrink fold changes for lowly expressed genes
# res <- lfcShrink(dds, contrast=contrast, res=res)

# # sort by p-value
# res <- res[order(res$padj),]
# # TODO explore IHW usage


# # store results
# svg(snakemake@output[["ma_plot"]])
# plotMA(res, ylim=c(-2,2))
# dev.off()

# write.table(as.data.frame(res), file=snakemake@output[["table"]])
