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

cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
cts <- cts[ , order(names(cts))]

coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE)
coldata <- coldata[order(row.names(coldata)), , drop=F]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=as.formula(snakemake@params[["model"]]))

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel=parallel)

# Write dds object as RDS
saveRDS(dds, file=snakemake@output[[1]])
# Write normalized counts
norm_counts = counts(dds, normalized=T)
write.table(data.frame("gene"=rownames(norm_counts), norm_counts), file=snakemake@output[[2]], sep='\t', row.names=F)
