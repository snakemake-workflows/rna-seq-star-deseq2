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
        adapters="-a {} ".format(config["trimming"]["adapter"]),
        others=config["params"]["cutadapt-pe"]
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    threads: 24
    wrapper:
        "0.66.0/bio/cutadapt/pe"


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
    threads: 24
    wrapper:
        "0.66.0/bio/cutadapt/se"
