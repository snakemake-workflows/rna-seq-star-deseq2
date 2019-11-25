def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        adapters_r1 = "-a {} {}".format(config["trimming"]["adapter"], config["params"]["cutadapt-pe"]),
        adapters_r2 = "-A {} {}".format(config["trimming"]["adapter"], config["params"]["cutadapt-pe"]),
        others = "--minimum-length 1"
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.36.0/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastq
    output:
        fastq="trimmed/{sample}-{unit}.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["trimming"]["adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.36.0/bio/cutadapt/se"
