rule count_matrix:
    input:
        expand(
            "star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample"].tolist(),
        strand=get_strandedness(units),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule deseq2_init:
    input:
        counts="counts/all.tsv",
    output:
        "deseq2/all.rds",
    params:
        samples=config["samples"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "deseq2/all.rds",
    output:
        report("results/pca.svg", "../report/pca.rst"),
    params:
        pca_labels=config["pca"]["labels"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log",
    script:
        "../scripts/plot-pca.R"


rule deseq2:
    input:
        "deseq2/all.rds",
    output:
        table=report(
            "results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"
        ),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast,
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log",
    threads: get_deseq2_threads
    script:
        "../scripts/deseq2.R"
