rule count_matrix:
    input:
        get_star_output,
    output:
        "results/counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule deseq2_init:
    input:
        counts="results/counts/all.tsv",
    output:
        "results/deseq2/all.rds",
    params:
        samples=config["samples"],
        model=config["diffexp"]["model"],
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log",
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"


rule pca:
    input:
        "results/deseq2/all.rds",
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
        "results/deseq2/all.rds",
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
