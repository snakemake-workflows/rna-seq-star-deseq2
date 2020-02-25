log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")


library("DESeq2")
library("pracma")
library("ashr")
library("IHW")
library("ggplot2")

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

IHWactive <- snakemake@params[["IHWactive"]]
IHWactive <- sapply(IHWactive,toupper)

type <- snakemake@params[["lfcShrinkType"]]
lfcShrink <- snakemake@params[["lfcShrink"]]
lfcShrink <- sapply(lfcShrink,toupper)

IHWalpha <- as.numeric(snakemake@params[["IHWalpha"]])
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
 		plotMA(res, alpha=alpha)
 		dev.off()

 		write.table(as.data.frame(res), file=snakemake@output[["table"]])

 		if (IHWactive) {
	 		res <- as.data.frame(res)
			# histogram 
			svg(snakemake@output[["pvalHisto1"]])
			g <- ggplot(res, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
			print(g)

			# IHW
			ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = IHWalpha)
			write.table(as.data.frame(ihwRes), file=snakemake@output[["IHWData"]])
			# estimatedWeights
			svg(snakemake@output[["IHWPlots"]])
			print(plot(ihwRes))
			dev.off()

			# decisionBoundary
			svg(snakemake@output[["IHWPlots2"]])
			print(plot(ihwRes, what = "decisionboundary"))
			dev.off()
		}
	}else{
		res <- results(dds, contrast=contrast, parallel=parallel, alpha=alpha)
		if (lfcShrink) {
 			res <- lfcShrink(dds, contrast=contrast, type=type,res=res)
 		}
 		# sort by p-value
 		res <- res[order(res$padj),]
 		# store results
 		svg(snakemake@output[["ma_plot"]])
 		plotMA(res, alpha=alpha)
 		dev.off()

 		write.table(as.data.frame(res), file=snakemake@output[["table"]])

 		if (IHWactive) {
	 		res <- as.data.frame(res)
			# histogram 
			svg(snakemake@output[["pvalHisto1"]])
			g <- ggplot(res, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
			print(g)

			# IHW
			ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = IHWalpha)
			write.table(as.data.frame(ihwRes), file=snakemake@output[["IHWData"]])
			# estimatedWeights
			svg(snakemake@output[["IHWPlots"]])
			print(plot(ihwRes))
			dev.off()

			# decisionBoundary
			svg(snakemake@output[["IHWPlots2"]])
			print(plot(ihwRes, what = "decisionboundary"))
			dev.off()
		}
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
 	plotMA(res, alpha=alpha)
 	dev.off()

 	write.table(as.data.frame(res), file=snakemake@output[["table"]])

 	if (IHWactive) {
	 	res <- as.data.frame(res)
		# histogram 
		svg(snakemake@output[["pvalHisto1"]])
		g <- ggplot(res, aes(x = pvalue)) + geom_histogram(binwidth = 0.025, boundary = 0)
		print(g)

		# IHW
		ihwRes <- ihw(pvalue ~ baseMean,  data = res, alpha = IHWalpha)
		write.table(as.data.frame(ihwRes), file=snakemake@output[["IHWData"]])
		# estimatedWeights
		svg(snakemake@output[["IHWPlots"]])
		print(plot(ihwRes))
		dev.off()

		# decisionBoundary
		svg(snakemake@output[["IHWPlots2"]])
		print(plot(ihwRes, what = "decisionboundary"))
		dev.off()
	}
}
graphics.off()