#TODO: Merge with deseq2.R or at least share common code

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

dds <- readRDS(snakemake@input[[1]])


#Parse interaction

environment <- snakemake@params[["interaction"]][["environment"]]
contrast <- snakemake@params[["interaction"]][["contrast"]]
cond <- contrast[1]
levelA <- contrast[2]
levelB <- contrast[3]

#Restructure the data to represent combinations

getColumn <- function (index,ds){
	return(ds[[index]])
}

#infer variables needed for model from provided environment and contrast

variables <- unlist(list(cond,names(environment))) #TODO: There is a more elegant solution, I'm certain

vars <- lapply(variables,getColumn,ds=dds)
#R ...
tmp <- do.call(paste0,vars)

newgroup <- factor(tmp)

dds$group <- newgroup
design(dds) <- ~ group
dds <- DESeq(dds)


#Construct levelA and levelB from environment and contrast

levelAcomb <- ""
levelBcomb <- ""


for (var in variables){
	if (exists(var, where=environment)){
		levelAcomb <- paste0(levelAcomb,environment[var])
		levelBcomb <- paste0(levelBcomb,environment[var])
	}
	else if (var == cond){
		levelAcomb <- paste0(levelAcomb,levelA)
		levelBcomb <- paste0(levelBcomb,levelB)
	}
	else{
		cat("Error:",var,": variable is not defined as either environment or contrast\n")
		stop("Terminating ...")
	}
}

cat("Level A:",levelAcomb,"\n")
cat("Level B:",levelBcomb,"\n")

contrast <- c("group",levelAcomb,levelBcomb)

cat("ResultNames:",resultsNames(dds),"\n")

res <- results(dds, contrast=contrast, parallel=parallel)

# shrink fold changes for lowly expressed genes #TODO Enable Shrinkage for Interactions?
#res <- lfcShrink(dds, contrast=contrast, res=res)

# sort by p-value
res <- res[order(res$padj),]
# TODO explore IHW usage


# store results
svg(snakemake@output[["ma_plot"]])
plotMA(res, ylim=c(-2,2))
dev.off()

write.table(as.data.frame(res), file=snakemake@output[["table"]])
