log <- file(snakemake@log[[1]], open = "wt")
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

counts_data <- read.table(
  snakemake@input[["counts"]],
  header = TRUE,
  row.names = "gene",
  check.names = FALSE
)
counts_data <- counts_data[, order(names(counts_data))]

col_data <- read.table(
  snakemake@params[["samples"]],
  header = TRUE,
  row.names = "sample_name",
  check.names = FALSE
)
col_data <- col_data[order(row.names(col_data)), , drop = FALSE]


coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample_name", check.names=FALSE)
coldata <- coldata[order(row.names(coldata)), , drop=F]

dds <- DESeqDataSetFromMatrix(countData=counts_data,
                              colData=col_data,
                              design=as.formula(snakemake@params[["model"]]))

# remove uninformative columns
dds <- dds[rowSums(counts(dds)) > 1, ]
# normalization and preprocessing
dds <- DESeq(dds, parallel = parallel)

# Write dds object as RDS
saveRDS(dds, file = snakemake@output[[1]])
# Write normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.table(
  data.frame(
    "gene" = rownames(norm_counts),
    norm_counts
  ),
  file = snakemake@output[[2]],
  sep = "\t",
  row.names = FALSE
)
