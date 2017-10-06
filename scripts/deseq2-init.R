library("DESeq2")

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
dds <- DESeqDataSetFromMatrix(countData=snakemake@input[["counts"]],
                              colData=snakemake@input[["samples"]],
                              design=~ condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]
# TODO optionally allow to collapse technical replicates
dds <- DESeq(dds)

save(dds, file=snakemake.output[[1]])
