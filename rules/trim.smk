def get_fastq(wildcards):
    columns = "fq1" if is_single_end(**wildcards) else ["fq1", "fq2"]
    return units.loc[(wildcards.sample, wildcards.unit), columns]


rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-pe"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastq
    output:
        fastq="trimmed/{sample}-{unit}.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/se"
