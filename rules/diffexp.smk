def get_strandness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]

rule count_matrix:
    input:
        expand("star/{unit.sample}-{unit.unit}/ReadsPerGene.out.tab", unit=units.itertuples())
    output:
        "counts/all.tsv"
    params:
        samples=units["sample"].tolist(),
        strand=get_strandness(units)
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

rule deseq2_init:
    input:
        counts="counts/all.tsv"
    output:
        dds="deseq2/all.rds",
        dds2="deseq2/all2.rds"
    params:
        samples=config["samples"],
        formula=config["diffexp"]["advanced"]["formula"]["design"],
        time=config["diffexp"]["advanced"]["time"],
        reduced=config["diffexp"]["advanced"]["reduced"],
        group=config["diffexp"]["advanced"]["formula"]["group"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/init.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2-init.R"

rule pca:
    input:
        "deseq2/all.rds"
    output:
        report("results/pca.svg", "../report/pca.rst")
    params:
        pca_labels=config["pca"]["labels"],
        transformation=config["pca"]["transformation"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/pca.log"
    script:
        "../scripts/plot-pca.R"

rule meanSdPlot:
    input:
        "deseq2/all.rds"
    output:
        report("results/meanSdPlot.svg", "../report/meanSdPlot.rst")
    params:
        transformation=config["meanSdPlot"]["transformation"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/meanSdPlot.log"
    script:
        "../scripts/meanSdPlot.R"

rule heatmap:
    input:
        "deseq2/all.rds"
    output:
        report("results/heatmap.svg", "../report/heatmap.rst")
    params:
        samples=config["samples"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/heatmap.log"
    script:
        "../scripts/heatmap.R"

rule dispersionPlot:
    input:
        "deseq2/all.rds"
    output:
        report("results/dispersionPlot.svg", "../report/dispersionPlot.rst")
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/dispersionPlot.log"
    script:
        "../scripts/dispersionPlot.R"

def get_contrast(wildcards):
    return config["diffexp"]["advanced"]["contrasts"][wildcards.contrast]

checkpoint deseq2:
    input:
        dds="deseq2/all.rds",
        dds2="deseq2/all2.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        samples=config["samples"],
        formula=config["diffexp"]["advanced"]["formula"]["design"],
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"
