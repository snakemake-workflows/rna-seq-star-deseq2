This workflow performs differential expression analysis on single- or paired-end RNA-seq data.
After adapter removal with `fastp <https://github.com/OpenGene/fastp>`_, reads were mapped and gene counts were generated with `STAR <https://github.com/alexdobin/STAR>`_.
Gene counts of replicated were summed up.
Integrated normalization and differential expression analysis was conducted with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ following standard procedure as outlined in the manual.
The differential expression results can be explored in interactive reports created with `Datavzrd <https://datavzrd.github.io>`_, which can be found under Differential expression.
