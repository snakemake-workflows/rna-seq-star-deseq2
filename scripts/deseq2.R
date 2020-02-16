log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("DESeq2")
library("pracma")
library("ashr")
library("IHW")

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

alpha <- as.numeric(snakemake@params[["alpha"]])
# create all contrasts if no formula is given
# otherwise use given formula and contrasts
if (strcmp(Reduce(paste, deparse(formula)),"~1")) { 
	# in case we have a complex contrast use dds2 with grouped formula
	# otherwise use simple contrasts
	if(strcmp(contrast[[1]],'group')){ 
		res <- results(dds2, contrast=contrast, parallel=parallel, alpha=alpha)
 		res <- lfcShrink(dds2, contrast=contrast, type="ashr", res=res)
 		# sort by p-value
 		res <- res[order(res$padj),]
 		# store results
 		svg(snakemake@output[["ma_plot"]])
 		plotMA(res, main="apeglm")
 		dev.off()
 		write.table(as.data.frame(res), file=snakemake@output[["table"]])
	}else{
		res <- results(dds, contrast=contrast, parallel=parallel, alpha=alpha)
 		res <- lfcShrink(dds, contrast=contrast, type="ashr",res=res)
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
 	res <- lfcShrink(dds, contrast=contrast, res=res) ## only if user wants it. Will throw "LFC shrinkage type='normal' not implemented for designs with interactions"
 	# sort by p-value
 	res <- res[order(res$padj),]
 	# store results
 	svg(snakemake@output[["ma_plot"]])
 	plotMA(res, main="apeglm")
 	dev.off()
 	write.table(as.data.frame(res), file=snakemake@output[["table"]])
}

write.table(as.data.frame(sum(res$padj < alpha)), file=snakemake@output[["sumPadj"]])
