log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("pracma")
library("DESeq2")

# load deseq2 data
dds <- readRDS(snakemake@input[[1]])


transformation <- snakemake@params[["transformation"]]

if(strcmp(transformation,"rlog")){
	rld <- rlog(dds, blind=FALSE)
	svg(snakemake@output[[1]])
	plotPCA(rld, intgroup=snakemake@params[["pca_labels"]])
}else{
	vsd <- vst(dds, blind=FALSE)
   	svg(snakemake@output[[1]])
	plotPCA(vsd, intgroup=snakemake@params[["pca_labels"]])
}


while (!is.null(dev.list()))  dev.off()