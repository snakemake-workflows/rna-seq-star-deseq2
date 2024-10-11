# Check if `snakemake` object exists, if not, define it manually
if (!exists("snakemake")) {
  snakemake <- list(
    input = list(counts = "results/counts/all.tsv"),
    output = list("results/deseq2/all.rds", "results/deseq2/normcounts.tsv"),
    log = list("logs/deseq2/init.log"),  # Correct the log definition
    threads = 1,  # Or set the number of threads you want to use
    config = list(
      samples = "path/to/samples.tsv",
      diffexp = list(
        variables_of_interest = list(
          "some_variable" = list(base_level = "some_level")
        ),
        batch_effects = c("effect1", "effect2")
      )
    )
  )
}


# Redirecting logs
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(stringr)
library("DESeq2")

# Setup parallelization
parallel <- FALSE
if (snakemake@threads > 1) {
  library("BiocParallel")
  register(MulticoreParam(snakemake@threads))
  parallel <- TRUE
}

# Print the input file names for troubleshooting
cat("Counts file: ", snakemake@input[["counts"]], "\n")
cat("Samples file: ", snakemake@config[["samples"]], "\n")

# Read input files
counts_data <- read.table(
  snakemake@input[["counts"]],
  header = TRUE,
  row.names = "gene",
  check.names = FALSE
)
counts_data <- counts_data[, order(names(counts_data))]

col_data <- read.table(
  snakemake@config[["samples"]],
  header = TRUE,
  row.names = "sample_name",
  check.names = FALSE
)
col_data <- col_data[order(row.names(col_data)), , drop = FALSE]

# Set base level for variables of interest
for (vof in names(snakemake@config[["diffexp"]][["variables_of_interest"]])) {
  base_level <- snakemake@config[["diffexp"]][["variables_of_interest"]][[vof]][["base_level"]]
  col_data[[vof]] <- relevel(factor(col_data[[vof]]), base_level)
}

# Handle batch effects
batch_effects <- snakemake@config[["diffexp"]][["batch_effects"]]
for (effect in batch_effects) {
  if (str_length(effect) > 0) {
    col_data[[effect]] <- factor(col_data[[effect]])
  }
}



# build up formula with additive batch_effects and all interactions between the
# variables_of_interes

design_formula <- snakemake@config[["diffexp"]][["model"]]

if (str_length(design_formula) == 0) {
  batch_effects <- str_flatten(batch_effects, " + ")
  if (str_length(batch_effects) > 0) {
    batch_effects <- str_c(batch_effects, " + ")
  }
  vof_interactions <- str_flatten(
    names(snakemake@config[["diffexp"]][["variables_of_interest"]]),
    " * "
  )
  design_formula <- str_c("~", batch_effects, vof_interactions)
}

print(paste("Dimensions of counts_data: ", dim(counts_data)))
print(paste("Dimensions of col_data: ", dim(col_data)))

dds <- DESeqDataSetFromMatrix(
  countData = counts_data,
  colData = col_data,
  design = as.formula(design_formula)
)

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
