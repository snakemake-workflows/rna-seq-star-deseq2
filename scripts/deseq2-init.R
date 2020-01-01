log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("tidyverse")

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

# get the design formula specified in config
formula <- as.formula(snakemake@params[["formula"]])

# create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData=cts,
                              colData=coldata,
                              design=formula)

# group?
group <- snakemake@params[["group"]]
group <- apply(group,toupper)

if (group) {
	variables <- labels(terms(formula)) %>%
  					 strsplit('[:*]') %>%
  					 unlist()
  	#^ might be useful cause it returns a string with the columns, but problematic is to automate something like dds$...
  	# 2nd option shows how to do it with contrasts, would need to think about a possible implementation

 	#dds$group <- factor(paste0(dds$condition,dds$time))
	# design(dds) <- ~ group

	# dds$group <- factor(paste0(c("treated","untreated","treated","untreated"),c("4h","4h","8h","8h")))
	# design(dds) <- ~ group

	#dds <- DESeq(dds)
}




# remove uninformative rows with no counts 
dds <- dds[ rowSums(counts(dds)) > 1, ]

# is this a Time Course Experiment?
time <- snakemake@params[["time"]]
time <- apply(time,toupper)

if (time) {
    # Time Course Experiment
	reduced <- as.formula(snakemake@params[["reduced"]])
    ddsTCE_LRT <- DESeq(dds, test="LRT", reduced = reduced, parallel=parallel)

    # create in deseq2.R
    #resTCE_LRT <- results(ddsTCE_LRT)
}else{
   	dds <- DESeq(dds, parallel=parallel)
}



#new but might be unnecessary, would have to do it before the DESeq call
#try standard: parametric fitType
# suppressWarnings(tryDESeq <- try(DESeq(dds, parallel=parallel, quiet=TRUE), silent=TRUE))

# if(inherits(tryDESeq,"try-error")){
# 	# try local fitType
# 	suppressWarnings(tryDESeq <- try(DESeq(dds, fitType="local", parallel=parallel, quiet=TRUE), silent=TRUE))


# 	if(inherits(tryDESeq,"try-error")){
# 		# try mean fitType
# 		suppressWarnings(tryDESeq <- try(DESeq(dds, fitType="mean", parallel=parallel, quiet=TRUE), silent=TRUE))
# 	}

# 	if(inherits(tryDESeq,"try-error")){
# 	#all gene-wise dispersion estimates are within 2 orders of magnitude from the minimum value
# 	#so we will use gene-wise estimates as final estimates
# 		dds <- estimateSizeFactors(dds)
# 		dds <- estimateDispersionsGeneEst(dds)
# 		dispersions(dds) <- mcols(dds)$dispGeneEst
# 		tryDESeq <- nbinomWaldTest(dds)
# 		#try previous: try(nbinomWaldTest(dds),silent=TRUE)
# 	}
# }

saveRDS(tryDESeq, file=snakemake@output[[1]])



# should save the created coefs into deseq2/coefs.txt
coefs <- resultsNames(tryDESeq)
coefsList <- as.list(coefs)

coefsF <- file(snakemake@output[[2]], open="wt")
sink(coefsF)
 
for (i in 2:length(list)) { 
	print(coefsList[[i]])
	print("\n")
}
sink()