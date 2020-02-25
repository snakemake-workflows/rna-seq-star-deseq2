log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("pracma")
library("vsn")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])


transformation <- snakemake@params[["transformation"]]

if(strcmp(transformation,"norm")){
	# this gives log2(n + 1)
	ntd <- normTransform(dds)
   	svg(snakemake@output[[1]])
	meanSdPlot(assay(ntd))
	dev.off()
}

if(strcmp(transformation,"rlog")){
	rld <- rlog(dds)
	svg(snakemake@output[[1]])
	meanSdPlot(assay(rld))
	dev.off()
}

if(strcmp(transformation,"VST")){
	vsd <- vst(dds)
   	svg(snakemake@output[[1]])
	meanSdPlot(assay(vsd))
	dev.off()
}


