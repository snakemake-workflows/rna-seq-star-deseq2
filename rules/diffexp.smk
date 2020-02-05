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
        "deseq2/all.rds"
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

def get_contrast(wildcards):
    return config["diffexp"]["advanced"]["contrasts"][wildcards.contrast]

checkpoint deseq2_ALL:
    input:
        "deseq2/all.rds"
    output:
        directory("results/diffexpAll/")
    params:
        samples=config["samples"],
        ALL=config["diffexp"]["standard"]["createContrasts"]
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/deseq2_ALL.diffexp.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2_ALL.R"

checkpoint deseq2:
    input:
        "deseq2/all.rds"
    output:
        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
    params:
        contrast=get_contrast
    conda:
        "../envs/deseq2.yaml"
    log:
        "logs/deseq2/{contrast}.diffexp.log"
    threads: get_deseq2_threads()
    script:
        "../scripts/deseq2.R"

def aggregate_input(wildcards):
    checkpoint_output = checkpoints.deseq2_ALL.get(**wildcards).output[0]
    return expand(["results/diffexpAll/{i}.diffexp.tsv",
            "results/diffexpAll/{i}.ma-plot.svg"],
           i=glob_wildcards(os.path.join(checkpoint_output, '{i}.ma-plot.svg')).i)


rule addALLtoReport:
    input:
        tsv1="results/diffexpAll/{i}.diffexp.tsv",
        svg1="results/diffexpAll/{i}.ma-plot.svg"
    output:
        tsv2=report("results/diffexpAll2/{i}.diffexp.tsv", "../report/diffexpAll.rst"),
        svg2=report("results/diffexpAll2/{i}.ma-plot.svg", "../report/maAll.rst")
    params:
        i=lambda w: aggregate_input(w)
    run:
        shell("cp {input.tsv1} {output.tsv2} ")
        shell("cp {input.svg1} {output.svg2} ")
