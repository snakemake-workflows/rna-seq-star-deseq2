rule get_sra:
    output:
        "sra/{accession}_1.fastq",
        "sra/{accession}_2.fastq",
    log:
        "logs/get-sra/{accession}.log",
    wrapper:
        "v3.5.3/bio/sra-tools/fasterq-dump"


rule fastp_se:
    input:
        sample=get_units_fastqs,
    output:
        trimmed="results/trimmed/{sample}/{sample}-{unit}_single.fastq.gz",
        failed="results/trimmed/{sample}/{sample}-{unit}_single.failed.fastq.gz",
        html=report(
            "results/trimmed/{sample}/{sample}-{unit}_single.html",
            caption="../report/fastp.rst",
            category="quality control",
            subcategory="fastp",
            labels={
                "sample-unit": "{sample}-{unit}",
            },
        ),
        json="results/trimmed/{sample}/{sample}-{unit}.json",
    log:
        "logs/trimmed/{sample}/{sample}-{unit}.log",
    params:
        adapters=lookup(
            within=units,
            query="sample_name == '{sample}' & unit_name == '{unit}'",
            cols="fastp_adapters",
            default="",
        ),
        extra=lookup(
            within=units,
            query="sample_name == '{sample}' & unit_name == '{unit}'",
            cols="fastp_extra",
            default="--trim_poly_x --poly_x_min_len 7 --trim_poly_g --poly_g_min_len 7",
        ),
    threads: 4
    wrapper:
        "v7.1.0/bio/fastp"


rule fastp_pe:
    input:
        sample=get_units_fastqs,
    output:
        trimmed=[
            "results/trimmed/{sample}/{sample}-{unit}_R1.fastq.gz",
            "results/trimmed/{sample}/{sample}-{unit}_R2.fastq.gz",
        ],
        # Unpaired reads separately
        unpaired1="results/trimmed/{sample}/{sample}-{unit}.unpaired.R1.fastq.gz",
        unpaired2="results/trimmed/{sample}/{sample}-{unit}.unpaired.R2.fastq.gz",
        failed="results/trimmed/{sample}/{sample}-{unit}_paired.failed.fastq.gz",
        html=report(
            "results/trimmed/{sample}/{sample}-{unit}.html",
            caption="../report/fastp.rst",
            category="quality control",
            subcategory="fastp",
            labels={
                "sample-unit": "{sample}-{unit}",
            },
        ),
        json="results/trimmed/{sample}/{sample}-{unit}.json",
    log:
        "logs/trimmed/{sample}/{sample}-{unit}.log",
    params:
        adapters=lookup(
            within=units,
            query="sample_name == '{sample}' & unit_name == '{unit}'",
            cols="fastp_adapters",
            default="--detect_adapter_for_pe",
        ),
        extra=lookup(
            within=units,
            query="sample_name == '{sample}' & unit_name == '{unit}'",
            cols="fastp_extra",
            default="--trim_poly_x --poly_x_min_len 7 --trim_poly_g --poly_g_min_len 7",
        ),
    threads: 8
    wrapper:
        "v7.1.0/bio/fastp"
