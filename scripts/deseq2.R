log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("DESeq2")
library("pracma")
library("ashr")
library("apeglm")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}


coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

dds <- readRDS(snakemake@input[["dds"]])
dds2 <- readRDS(snakemake@input[["dds2"]])
# get contrasts
contrast <- c(snakemake@params[["contrast"]])

formula <- as.formula(snakemake@params[["formula"]])

type <- snakemake@params[["lfcShrinkType"]]
lfcShrink <- snakemake@params[["lfcShrink"]]
lfcShrink <- sapply(lfcShrink,toupper)

alpha <- as.numeric(snakemake@params[["alpha"]])
# create all contrasts if no formula is given
# otherwise use given formula and contrasts
if (strcmp(Reduce(paste, deparse(formula)),"~1")) { 
	# in case we have a complex contrast use dds2 with grouped formula
	# otherwise use simple contrasts
	if(strcmp(contrast[[1]],'group')){ 
		res <- results(dds2, contrast=contrast, parallel=parallel, alpha=alpha)

		if (lfcShrink) {
 			res <- lfcShrink(dds2, contrast=contrast, type=type, res=res)
 		}
 		# sort by p-value
 		res <- res[order(res$padj),]
 		# store results
 		svg(snakemake@output[["ma_plot"]])
 		plotMA(res, main="apeglm")
 		dev.off()
 		write.table(as.data.frame(res), file=snakemake@output[["table"]])
	}else{
		res <- results(dds, contrast=contrast, parallel=parallel, alpha=alpha)
		if (lfcShrink) {
 			res <- lfcShrink(dds, contrast=contrast, type=type,res=res)
 		}
 		# sort by p-value
 		res <- res[order(res$padj),]
 		# store results
 		svg(snakemake@output[["ma_plot"]])
 		plotMA(res, main="apeglm")
 		dev.off()
 		write.table(as.data.frame(res), file=snakemake@output[["table"]])
	}
}else{	

 	res <- results(dds, contrast=contrast, parallel=parallel, alpha=alpha)
 	if (lfcShrink) {
 		res <- lfcShrink(dds, contrast=contrast, res=res, type=type)
 	}
 	# sort by p-value
 	res <- res[order(res$padj),]
 	# store results
 	svg(snakemake@output[["ma_plot"]])
 	plotMA(res, main="apeglm")
 	dev.off()
 	write.table(as.data.frame(res), file=snakemake@output[["table"]])
}
