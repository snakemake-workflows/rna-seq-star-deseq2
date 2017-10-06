rule count_matrix:
    input:
        expand("star/{sample}/ReadsPerGene.out.tab", sample=samples.index)
    output:
        "counts/all.tsv"
    params:
        samples=samples.index
    script:
        "../scripts/count-matrix.py"


rule deseq2_init:
    input:
        counts="counts/all.tsv",
        samples="samples.tsv"
    output:
        "deseq2/all.RData"
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "deseq2/all.RData"
    output:
        "results/pca.pdf"
    params:
        pca_labels=config["pca"]["labels"]
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/pca.R"


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]


rule deseq2:
    input:
        "deseq2/all.RData"
    output:
        table="results/diffexp/{contrast}.diffexp.tsv",
        ma_plot="results/diffexp/{contrast}.ma-plot.pdf",
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    script:
        "../scripts/deseq2.R"
