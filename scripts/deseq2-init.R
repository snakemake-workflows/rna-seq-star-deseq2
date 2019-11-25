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

# colData and countData must have the same sample order, but this is ensured
# by the way we create the count matrix
cts <- read.table(snakemake@input[["counts"]], header=TRUE, row.names="gene", check.names=FALSE)
coldata <- read.table(snakemake@params[["samples"]], header=TRUE, row.names="sample", check.names=FALSE)

#formula <- snakemake@params[["contrast"]]

dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=~ type + condition)

# remove uninformative columns
dds <- dds[ rowSums(counts(dds)) > 1, ]


# normalization and preprocessing
#old
#dds <- DESeq(dds, parallel=parallel)

#new
#try standard: parametric fitType
suppressWarnings(tryDESeq <- try(DESeq(dds, parallel=parallel, quiet=TRUE), silent=TRUE))

if(inherits(tryDESeq,"try-error")){
	# try local fitType
	suppressWarnings(tryDESeq <- try(DESeq(dds, fitType="local", parallel=parallel, quiet=TRUE), silent=TRUE))


	if(inherits(tryDESeq,"try-error")){
		# try mean fitType
		suppressWarnings(tryDESeq <- try(DESeq(dds, fitType="mean", parallel=parallel, quiet=TRUE), silent=TRUE))
	}

	if(inherits(tryDESeq,"try-error")){
	#all gene-wise dispersion estimates are within 2 orders of magnitude from the minimum value
	#so we will use gene-wise estimates as final estimates
		dds <- estimateSizeFactors(dds)
		dds <- estimateDispersionsGeneEst(dds)
		dispersions(dds) <- mcols(dds)$dispGeneEst
		tryDESeq <- nbinomWaldTest(dds)
		#try previous: try(nbinomWaldTest(dds),silent=TRUE)
	}
}

saveRDS(tryDESeq, file=snakemake@output[[1]])
