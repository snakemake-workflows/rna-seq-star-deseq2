import glob

import pandas as pd
from snakemake.remote import FTP
from snakemake.utils import validate

ftp = FTP.RemoteProvider()

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str}).set_index("sample_name", drop=False).sort_index()

def get_final_output():
    final_output = expand("results/diffexp/{contrast}.diffexp.tsv",
                        contrast=config["diffexp"]["contrasts"])
    return final_output

validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str}).set_index(["sample_name", "unit_name"], drop=False).sort_index()
validate(units, schema="../schemas/units.schema.yaml")

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
    if not is_activated("mergeReads"):
        if config["trimming"]["activate"]:
            return expand("results/trimmed/{sample}_{unit}_R1.fastq.gz", 
                unit=units.loc[wildcards.sample, "unit_name"], sample=wildcards.sample)
        unit = units.loc[wildcards.sample]
        if all(pd.isna(unit["fq1"])):
            # SRA sample (always paired-end for now)
            accession = unit["sra"]
            return expand("sra/{accession}_R1.fastq", accession=accession)
        sample_units = units.loc[wildcards.sample]
        return sample_units["fq1"]
    if is_paired_end(wildcards.sample):
        return "results/merged/{sample}_R1.fastq.gz"
    return "results/merged/{sample}_single.fastq.gz"

def get_map_reads_input_R2(wildcards):
    if is_paired_end(wildcards.sample):
        if not is_activated("mergeReads"):
            if config["trimming"]["activate"]:
                return expand("results/trimmed/{sample}_{unit}_R1.fastq.gz", 
                    unit=units.loc[wildcards.sample, "unit_name"], sample=wildcards.sample)
            unit = units.loc[wildcards.sample]
            if all(pd.isna(unit["fq1"])):
                # SRA sample (always paired-end for now)
                accession = unit["sra"]
                return expand("sra/{accession}_R2.fastq", accession=accession)
            sample_units = units.loc[wildcards.sample]
            return sample_units["fq2"]
        return "results/merged/{sample}_R2.fastq.gz",
    return ""

def get_star_output(wildcards):
    res = []
    for unit in units.itertuples():
        if is_paired_end(unit.sample_name):
            lib = "pe"
        else:
            lib = "se"
        print(lib)
        res.append("results/star/{}/{}-{}/ReadsPerGene.out.tab".format(
            lib, unit.sample_name, unit.unit_name
        ))
    return res

def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list=["none"]
        return strand_list*units.shape[0]

def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6

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

def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]

def get_arriba_input(wc):
    if is_paired_end(wc.sample):
        return("results/star/{}/{}-{}/Aligned.out.bam".format(
            lib, wc.sample, wc.unit
        ))