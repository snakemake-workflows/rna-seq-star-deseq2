import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str}).set_index("sample_name", drop=False).sort_index()

def get_final_output():
    if config["diffexp"]["activate"]:
        final_output = expand("results/diffexp/{contrast}.diffexp.tsv",
                        contrast=config["diffexp"]["contrasts"])
    else:
        final_output = "results/counts/all.tsv",
    return final_output

validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str}).set_index(["sample_name", "unit_name"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

def get_fusion_output():
    if config["fusion"]["arriba"]["activate"]:
        fusion_output = expand("results/fusion/arriba/{sample}.fusions.tsv", sample=samples.sample_name)
    else:
        fusion_output = ""
    return fusion_output

def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(S=unit.sample_name, U=unit.unit_name, E=ending)
    else:
        # paired end local sample
        return expand("pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(S=unit.sample_name, U=unit.unit_name, E=ending), read=["fq1","fq2"])


def get_cutadapt_pipe_input(wildcards):
    files = list(sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq])))
    assert(len(files) > 0)
    return files


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert all_single or all_paired, "invalid units for sample {}, must be all paired end or all single end".format(sample)
    return all_paired

def get_map_reads_input_R1(wildcards):
    if is_paired_end(wildcards.sample):
        return "results/merged/{sample}_R1.fastq.gz"
    return "results/merged/{sample}_single.fastq.gz"

def get_map_reads_input_R2(wildcards):
    if is_paired_end(wildcards.sample):
        return "results/merged/{sample}_R2.fastq.gz",
    return ""

def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))

def get_fastqs(wc):
    if config["trimming"]["activate"]:
        return expand("results/trimmed/{sample}/{unit}_{read}.fastq.gz", 
            unit=units.loc[wc.sample, "unit_name"], sample=wc.sample, read=wc.read)
    unit = units.loc[wc.sample]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=wc.read[-1])
    fq = "fq{}".format(wc.read[-1]) 
    return units.loc[wc.sample, fq].tolist()