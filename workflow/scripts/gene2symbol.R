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

annotated <- merge(df, g2g, by.x="gene", by.y="ensembl_gene_id")
annotated$gene <- ifelse(annotated$external_gene_name == '', annotated$gene, annotated$external_gene_name)
annotated$external_gene_name <- NULL
write.table(annotated, snakemake@output[["symbol"]], sep='\t', row.names=F)
