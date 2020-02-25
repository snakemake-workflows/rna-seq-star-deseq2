# Guidance for the config.yaml
Written by Julian Kremer (@jukre111)

# This file was made to guide you through the config.yaml. The order is: 1. **Set up** 2. **Plot Options** 3. **Differential Expression Options**

## Set up

##### samples: Path to sample sheet
##### units: Path to units sheet

#### trimming
Options:
- skip: disable or enable trimming
    - true
    - false


- adapter: sequencing adapter, e.g. ACGGATCGATCGATCGATCGAT


#### ref
Options:
- index: Path to STAR index


- annotation: Path to gtf-File

- Detail:
    - Possible command to generate the index if not given:
    - STAR --runMode genomeGenerate --runThreadN "NumberOfCores" --genomeDir "PathToSTARIndex" --genomeFastaFiles "pathToFastaFiles"
    - Problems occuring in this step are most likely related to STAR

For more detail: [STAR-Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) and [STAR-Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)


## Plot Options

### PCA
Options:
- activate: enable or disable PCA
    - true
    - false


- transformation: transformation on the count data
    - VST
    - rlog 


- labels: columns of the sample sheet to use for PCA
    - column
- Detail:
    - VST: Variance-stabilizing transformation (needs at least 1000 genes)
    - rlog: Regularized log transformation
    
For more detail: [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

### meanSdPlot
Options:
- activate: enable or disable meanSdPlot
    - true
    - false


- transformation: transformation on the count data
    - VST
    - rlog 
    - norm

- Detail:
    - VST: Variance-stabilizing transformation (needs at least 1000 genes)
    - rlog: Regularized log transformation
    - norm: Log2 with a pseudocount of 1
    
For more detail: [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) and [vsn-Package](https://bioconductor.org/packages/release/bioc/html/vsn.html)


### heatmap
Option:
- activate: enable or disable heatmap
    - true
    - false

- Detail:
    - Heatmap of the count matrix

For more detail: [DESeq2-Vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-quality-assessment-by-sample-clustering-and-visualization) and [pheatmap-Package](https://cran.r-project.org/web/packages/pheatmap/index.html)


### dispersionPlot
Option:
- activate: enable or disable dispersionPlot
    - true
    - false

- Detail:
    - Final estimates shrunk from the gene-wise estimates towards the fitted estimates

For more detail: [DESeq2-Vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#data-quality-assessment-by-sample-clustering-and-visualization) and [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)


### IndependentHypothesisWeighting
Options:
- activate: enable or disable IndependentHypothesisWeighting
    - true
    - false


- alpha: nominal level for FDR control
    - numeric value, e.g. 0.1


- Detail:
    - Multiple testing procedure. Takes p-values, covariates, significance level and then calculates weights for each p-value

For more detail: [IHW-Package](http://bioconductor.org/packages/release/bioc/html/IHW.html) and [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)


## Differential Expression Options

### Standard:  createContrasts (for inexperienced users)
Option:
- createContrasts: enable or disable createContrasts
    - true
    - false


- Detail:
    - Automatically create contrasts for inexperienced users. n(n-1)/2 combinations for each column, where n is the given number of different effects in a column. E.g. column **"condition"** has the effects **"control"**,**"mutant1"** and **"mutant2"** (n=3) and thus will create the contrasts **"control-vs-mutant1"**, **"control-vs-mutant2"** and **"mutant1-vs-mutant2"**
    - Additionally create a linear combination of effects among all columns in a new notional column called  **"group"**. The different contrasts will also be determined for the **"group"** column.

For more detail: [DESeq2-Package](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) and for the creation, estimateContrasts() in the Snakefile

### Advanced: (for experienced users)

#### lfcShrink 
Options:
- activate: enable or disable lfcShrink
    - true
    - false


- type: Log2 fold change shrinkage for visualization and ranking of genes
    - normal
    - ashr


- Detail:
    - normal: Standard log2 fold change shrinkage, design is **not** allowed to have interaction terms (: or * in the design)
    - ashr: MLE (maximum likelihood estimates) log2 fold changes and standard errors from DESeq used by ashr-Package to create Empirical Bayes approach for large-scale hypothesis testing and false discovery rate (FDR) estimation

For more detail: [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) and [ashr-Package](https://cran.r-project.org/web/packages/ashr/index.html)


#### alpha 
Option:
- alpha: always enabled
    - numeric value, e.g. 0.1


- Detail:
    - The value will be used in DESeq2 results(...)-function and MAplot

For more detail: [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)

#### time 
Options:
- activate: enable or disable Time Course Experiment
    - true
    - false


- reduced: reduced model for LRT (likelihood ratio test)
    - design formula, e.g. ~ condition*condition2
    
    
- threshold: Threshold for the color key range of the given heatmap
    - integer, e.g. 5

- Detail:
    - reduced: reduced design formula vs given design formula
    - E.g. ~condition+time+condition:time as given design formula and ~condition+time as reduced design formula results in a LRT with removed condition-specific differences over time. Genes with small p values from this test are those which at one or more time points after time 0 showed a condition-specific effect. Note therefore that this will not give small p values to genes that moved up or down over time in the same way in both condition effects.

For more detail: [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) and [rnaseqGene-Vignette](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments)


#### formula 
Options:
- design:
    - design formula, e.g. ~ condition*condition2
    
    
- group: enable or disable a grouped design
    - true
    - false



- Detail:
    - group: Create a linear combination of effects among all columns in a new notional column called  **"group"**. The different contrasts will have to be set under **contrasts**.

For more detail: [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) and [rnaseqGene-Vignette](https://www.bioconductor.org/packages/devel/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#time-course-experiments)


#### contrasts 
Options:
- arbitrary-name: 
    - column name
    - effect1
    - effect2


- Detail:
    - arbitrary-name: A name for the contrast
        - column name: Name for the column
        - effect1: effect of the column
        - effect2: another effect of the column
    - E.g.: 
        - Without **group**: 
         treated-vs-untreated:
            - condition
            - treated
            - untreated
        - With **group**: (let's assume we have another column **"time"**, with effects, **"4h"** and **"8h"**, respectively. _Column name can only be group_)
         4htreated-vs-4huntreated:
            - group
            - 4htreated
            - 4huntreated
        
For more detail: [DESeq2-Paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)
