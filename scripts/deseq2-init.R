log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("pracma")
library("DESeq2")
#library("tidyverse")

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

# get the design formula specified in config
formula <- as.formula(snakemake@params[["formula"]])

# create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=formula)

# remove uninformative rows with no counts 
dds <- dds[ rowSums(counts(dds)) > 1, ]

if (strcmp(Reduce(paste, deparse(formula)),"~1")) {
	formula <- do.call(paste, c(as.list(colnames(coldata)), sep = "*"))
	formula <- as.formula(paste0("~",formula))
	design(dds) <- formula
	dds <- DESeq(dds)

	dds2 <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=formula)
	dds2$group <- factor(do.call(paste0,coldata)) 										 	
	design(dds2) <- ~ group
	dds2 <- DESeq(dds2)
	saveRDS(dds2, file=snakemake@output[["dds2"]])
}

# group
group <- snakemake@params[["group"]]
group <- sapply(group,toupper)

if (group) {
	dds$group <- factor(do.call(paste0,coldata)) 										 	
	design(dds) <- ~ group
	dds <- DESeq(dds)
}

# is this a Time Course Experiment?
time <- snakemake@params[["time"]]
time <- sapply(time,toupper)

if (time) {
    # Time Course Experiment
	reduced <- as.formula(snakemake@params[["reduced"]])
    dds <- DESeq(dds, test="LRT", reduced = reduced, parallel=parallel)

    # create in deseq2.R
    #resTCE_LRT <- results(ddsTCE_LRT)
}else{
   	dds <- DESeq(dds, parallel=parallel)
}

saveRDS(dds, file=snakemake@output[["dds"]])