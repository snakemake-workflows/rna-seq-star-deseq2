log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("apeglm")
library("pracma")
library("gtools")
library("dplyr")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}


coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

dds <- readRDS(snakemake@input[[1]])

ALL <- snakemake@params[["ALL"]]
ALL <- sapply(ALL,toupper)

if (ALL){

	# create most useful non group contrasts
	for (i in 1:length(colnames(coldata))) {
		contrast <- c(colnames(coldata)[i],coldata[i] %>% sapply(levels))
    	res <- results(dds, contrast=contrast, parallel=parallel)
    	svg(paste0(snakemake@output[[1]], do.call(paste0, c(as.list(contrast))),".ma-plot.svg"))
    	plotMA(res, ylim=c(-2,2))
    	dev.off()
		write.table(as.data.frame(res), file=paste0(snakemake@output[[1]], do.call(paste0, c(as.list(contrast))),".diffexp.tsv"))
  	}


  	# create most useful group contrasts
  	dds$group <- factor(do.call(paste0,coldata))	
	design(dds) <- ~ group
	dds <- DESeq(dds)

	# all possible combinations of levels
	contrasts <- combinations(n=length(levels(dds$group)),r=ncol(coldata),c(levels(dds$group)),repeats.allowed = FALSE)
	for (i in 1:nrow(contrasts)) {
		contrast <- c("group",contrasts[i,1:ncol(coldata)])
    	res <- results(dds, contrast=contrast, parallel=parallel)
    	svg(paste0(snakemake@output[[1]], do.call(paste0, c(as.list(contrast))),".ma-plot.svg"))
    	plotMA(res, ylim=c(-2,2))
    	dev.off()
		write.table(as.data.frame(res), file=paste0(snakemake@output[[1]], do.call(paste0, c(as.list(contrast))),".diffexp.tsv"))
   	}
}