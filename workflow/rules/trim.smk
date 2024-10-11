rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    threads: 32
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "v3.5.3/bio/sra-tools/fasterq-dump"


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/cutadapt/{sample}/{unit}.{fq}.{ext}"),
    log:
        "logs/pipe-fastqs/catadapt/{sample}_{unit}.{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fastq|fastq\.gz",
    threads: 1  ## this does something special when running using pipe() output directives
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="results/trimmed/{sample}_{unit}_R1.fastq.gz",
        fastq2="results/trimmed/{sample}_{unit}_R2.fastq.gz",
        qc="results/trimmed/{sample}_{unit}.paired.qc.txt",
    log:
        "logs/cutadapt/{sample}_{unit}.log",
    params:
        extra=config["params"]["cutadapt-pe"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 64
    wrapper:
        "v3.5.3/bio/cutadapt/pe"


rule cutadapt_se:
    input:
        get_cutadapt_input,
    output:
        fastq="results/trimmed/{sample}_{unit}_single.fastq.gz",
        qc="results/trimmed/{sample}_{unit}_single.qc.txt",
    log:
        "logs/cutadapt/{sample}_{unit}.log",
    params:
        extra=config["params"]["cutadapt-se"],
        adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 64
    wrapper:
        "v3.5.3/bio/cutadapt/se"
