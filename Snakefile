import pandas as pd
from snakemake.utils import validate, min_version
import itertools
##### set minimum snakemake version #####
min_version("5.1.2")
shell.executable("/bin/bash")

##### load config and sample sheets #####

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"]).set_index("sample", drop=False).astype(str)
validate(samples, schema="schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample", "unit"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="schemas/units.schema.yaml")

simpleContrasts={}
complexContrasts={}

def getSimpleContrasts(columns,rows):
    simpleContrasts={}
    i=0
    # simple Contrasts e.g. treated vs untreated
    for key in columns: # n*(n-1)/2 combinations
        for values in itertools.combinations(rows[i],2):
            if values[0]!=values[1]:
                name=values[0]+"-vs-"+values[1]
                simpleContrasts[name]=[key,values[0],values[1]]
        i+=1
    return simpleContrasts

def getComplexContrasts(columns):
    complexRows=[]
    complexContrasts={}
    for index,row in samples.iterrows():
        string=""
        for col in range(0,len(columns)): # concatenate rows
            string+=row[1::][col]
        complexRows.append(string)
    # complex Contrasts e.g. 4hourstreated vs 8hourstreated
    for values in itertools.combinations(complexRows,2): # n*(n-1)/2 combinations
        if values[0]!=values[1]:
            name=values[0]+"-vs-"+values[1]
            complexContrasts[name]=['group',values[0],values[1]]
    return complexContrasts

def estimateContrasts():
    columns=list(samples.keys())[1::] # get column names without overhead
    rows=[]
    for key in columns:
        rows.append(list(samples.get(key).unique())) # get distinct rows values

    formula=config["diffexp"]["advanced"]["formula"]["design"]
    if formula=="~ 1":
        global simpleContrasts
        global complexContrasts
        simpleContrasts=getSimpleContrasts(columns,rows)
        complexContrasts=getComplexContrasts(columns)
        print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        print("::: No custom formula design was choosen, trying to create all useful contrasts :::")
        print(":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::")
        config["diffexp"]["advanced"]["contrasts"]={**simpleContrasts,**complexContrasts}
    else:
        print("::::::::::::::::::::::::::::::::::::::::::::::")
        print("::: Your custom formula design was choosen :::")
        print("::::::::::::::::::::::::::::::::::::::::::::::")

def additionalInput():
    input=[]
    if config["meanSdPlot"]["activate"]:
        input.extend(expand("results/meanSdPlot.svg"))
    if config["heatmap"]["activate"]:
        input.extend(expand("results/heatmap.svg"))
    if config["pca"]["activate"]:
        input.extend(expand("results/pca.svg"))
    if config["dispersionPlot"]["activate"]:
        input.extend(expand("results/dispersionPlot.svg"))
    if config["IndependentHypothesisWeighting"]["activate"]:
        input.extend(expand(["results/diffexp/IHW/{contrast}IHW.tsv",
                        "results/diffexp/IHW/{contrast}IHW.svg",
                        "results/diffexp/IHW/{contrast}IHW2.svg",
                        "results/diffexp/IHW/{contrast}pvalHisto.svg"],
                       contrast=config["diffexp"]["advanced"]["contrasts"]))
    if config["diffexp"]["advanced"]["time"]["activate"]:
        input.extend(expand("results/heatmapTime.svg"))
    return input

estimateContrasts()
##### target rules #####

rule all:
    input:
        expand(["results/diffexp/{contrast}.diffexp.tsv",
                "results/diffexp/{contrast}.ma-plot.svg"],
               contrast=config["diffexp"]["advanced"]["contrasts"]),
        additionalInput(),
        "qc/multiqc_report.html"


##### setup singularity #####

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
singularity: "docker://continuumio/miniconda3"


##### setup report #####

report: "report/workflow.rst"


##### load rules #####

include: "rules/common.smk"
include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/diffexp.smk"
include: "rules/qc.smk"
