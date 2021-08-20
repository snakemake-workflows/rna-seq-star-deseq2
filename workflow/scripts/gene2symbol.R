library(biomaRt)
library(tidyverse)

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because 
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while ( class(mart)[[1]] != "Mart" ) {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        stop(
        )
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "uswest",
                     uswest = "asia",
                     asia = "www",
                     www = {
                       # wait before starting another round through the mirrors,
                       # hoping that intermittent problems disappear
                       Sys.sleep(30)
                       "useast"
                     }
              )
    }
  )
}


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


