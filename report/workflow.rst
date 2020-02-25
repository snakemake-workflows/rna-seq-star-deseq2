This workflow performs differential expression analysis on single- or paired-end RNA-seq data.
After adapter removal with `Cutadapt <http://cutadapt.readthedocs.io>`_, reads were mapped and gene counts were generated with `STAR <https://github.com/alexdobin/STAR>`_.
Gene counts of replicated were summed up.
Integrated normalization and differential expression analysis was conducted with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_ following standard procedure as outlined in the manual.
