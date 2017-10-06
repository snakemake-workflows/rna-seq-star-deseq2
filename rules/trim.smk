def get_fastq(wildcards):
    return samples.loc[wildcards.sample, ["fq1", "fq2"]].dropna()


rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}.1.fastq.gz",
        fastq2="trimmed/{sample}.2.fastq.gz",
        qc="trimmed/{sample}.qc.txt"
    params:
        config["cutadapt"]["params"]
    log:
        "logs/cutadapt/{sample}.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastq
    output:
        fastq="trimmed/{sample}.fastq.gz",
        qc="trimmed/{sample}.qc.txt"
    params:
        config["cutadapt"]["params"]
    log:
        "logs/cutadapt/{sample}.log"
    wrapper:
        "0.17.4/bio/cutadapt/se"
