library("msigdbr")
library("GSVA")
library("pheatmap")

species_name <- snakemake@config[["ref"]][["species"]]
species_name <- gsub("(^|\\s)([a-z])", "\\1\\U\\2", tools::toTitleCase(gsub("_", " ", species_name)), perl=TRUE)

# Load just the Hallmark gene sets from MSigDB
gene_sets = msigdbr(species = species_name, category = 'H')

# Format and convert to GMT
gene_sets <- gene_sets[,c("gs_name","ensembl_gene")]
gene_sets <- gene_sets[!duplicated(gene_sets),]
gmt <- split(x = gene_sets$ensembl_gene, f = gene_sets$gs_name)

# Load the data
norm_counts <- read.table(snakemake@input[[1]], sep='\t', header=1, row.names = "gene")
norm_counts_log <- log2(norm_counts + 1)

# Run GSVA
gsva.param <- gsvaParam(as.matrix(norm_counts_log), gmt, kcdf="Poisson")
data.gsva <- gsva(gsva.param)

# Write it out; look at this table as a heatmap in python
write.table(data.gsva, 'hallmark_GSVA.csv', sep=',', col.names= NA, file = snakemake@output[[1]])

# Save heatmap 
svg(snakemake@output[[2]])
pheatmap(data.gsva,
         scale = "row",
         clustering_method = "complete", 
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",  
         show_rownames = TRUE,  
         show_colnames = TRUE, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50), 
         main = "GSVA Heatmap"
)
dev.off()
