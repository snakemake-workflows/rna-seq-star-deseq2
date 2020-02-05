log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("DESeq2")
#library("tidyverse")
library("apeglm")
library("pracma")
library("gtools")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

dds <- readRDS(snakemake@input[[1]])

# creates a plot for given coef
 contrast <- snakemake@params[["contrast"]]
# res <- results(dds, name="condition_untreated_vs_treated", parallel=parallel)
# res <- lfcShrink(dds, coef="condition_untreated_vs_treated", type="apeglm", res=res)
# res <- res[order(res$padj),]
# svg(snakemake@output[["ma_plot"]])
# plotMA(res, main="apeglm")
# dev.off()
# write.table(as.data.frame(res), file=snakemake@output[["table"]])



coefs <- resultsNames(dds)
coefsList <- as.list(coefs)

 
for (i in 1:length(coefsList)) {
	if(strcmp(coefsList[[i]],"Intercept")){

	}else{
		res <- results(dds, name=coefsList[[i]], parallel=parallel)
		res <- lfcShrink(dds, coef=coefsList[[i]], type="apeglm", res=res)
		res <- res[order(res$padj),]
		svg(snakemake@output[["ma_plot"]])
		plotMA(res, main="apeglm")
		dev.off()
		write.table(as.data.frame(res), file=snakemake@output[["table"]])
	}
}




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
