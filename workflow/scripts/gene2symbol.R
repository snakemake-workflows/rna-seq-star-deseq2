log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(biomaRt)

mart <- biomaRt::useEnsembl(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = paste0(snakemake@params[["species"]], "_gene_ensembl"),
  version = snakemake@params[["version"]]
)

df <- read.table(snakemake@input[["counts"]], sep='\t', header=1)

g2g <- biomaRt::getBM(
            attributes = c( "ensembl_gene_id",
                            "external_gene_name"),
            filters = "ensembl_gene_id",
            values = df$gene,
            mart = mart,
            )

df$symbol <- g2g$external_gene_name[match(df$gene, g2g$ensembl_gene_id)]
write.table(df, snakemake@output[["symbol"]], sep='\t', row.names=F)
