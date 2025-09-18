## RSEQC

FASTQC_WILDCARDS = []
for unit in units.itertuples():
    reads = ("R1", "R2") if is_paired_end(unit.sample_name) else ("R1",)
    for read in reads:
        FASTQC_WILDCARDS.append(
            {"sample": unit.sample_name, "unit": unit.unit_name, "read": read}
        )

FASTQC_ZIP_OUTPUTS = [
    "results/qc/fastqc/{sample}_{unit}_{read}_fastqc.zip".format(**entry)
    for entry in FASTQC_WILDCARDS
]

def get_fastqc_fastq(wildcards):
    fastqs = get_fq(wildcards)
    read_map = {"R1": "fq1", "R2": "fq2"}
    try:
        key = read_map[wildcards.read]
    except KeyError:
        raise ValueError(
            "Invalid read value '{read}' for sample {sample} unit {unit}".format(
                read=wildcards.read, sample=wildcards.sample, unit=wildcards.unit
            )
        )

    fastq = fastqs.get(key)
    if fastq is None:
        raise ValueError(
            "Read {read} not available for sample {sample} unit {unit}".format(
                read=wildcards.read, sample=wildcards.sample, unit=wildcards.unit
            )
        )

    return fastq

def get_multiqc_inputs(wildcards):
    inputs = list(FASTQC_ZIP_OUTPUTS)
    inputs.extend(
        expand(
            "results/star/{unit.sample_name}_{unit.unit_name}/Aligned.sortedByCoord.out.bam",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.junctionanno.junction.bed",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.junctionsat.junctionSaturation_plot.pdf",
            unit=units.itertuples(),
        )
    )
    inputs.extend(

        expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.infer_experiment.txt",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.stats.txt",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.inner_distance_freq.inner_distance.txt",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(

            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.readdistribution.txt",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.readdup.DupRate_plot.pdf",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(
            "results/qc/rseqc/{unit.sample_name}_{unit.unit_name}.readgc.GC_plot.pdf",
            unit=units.itertuples(),
        )
    )
    inputs.extend(
        expand(
            "logs/rseqc/rseqc_junction_annotation/{unit.sample_name}_{unit.unit_name}.log",
            unit=units.itertuples(),
        )
    )
    return inputs


rule fastqc:
    input:
        fastq=get_fastqc_fastq,
    output:
        html="results/qc/fastqc/{sample}_{unit}_{read}_fastqc.html",
        zip="results/qc/fastqc/{sample}_{unit}_{read}_fastqc.zip",
    threads: 4
    log:
        "logs/fastqc/{sample}_{unit}_{read}.log",
    params:
        extra="",
    wrapper:
        "v3.5.3/bio/fastqc"


rule rseqc_gtf2bed:
    input:
        "resources/genome.gtf",
    output:
        bed="results/qc/rseqc/annotation.bed",
        db=temp("results/qc/rseqc/annotation.db"),
    threads: 32
    log:
        "logs/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionanno.junction.bed",
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log[0]} 2>&1;
        """


rule rseqc_junction_saturation:
    input:
        bam="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.junctionsat.junctionSaturation_plot.pdf",
    priority: 1
    threads: 32
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}_{unit}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        """junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1;
	"""


rule rseqc_stat:
    input:
        "results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.stats.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}_{unit}.log",
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        bam_stat.py -i {input} > {output} 2> {log};
        """


rule rseqc_infer:
    input:
        bam="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.infer_experiment.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}_{unit}.log",
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log};
        """


rule rseqc_innerdis:
    input:
        bam="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    output:
        "results/qc/rseqc/{sample}_{unit}.inner_distance_freq.inner_distance.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    threads: 32
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1;
        """


rule rseqc_readdis:
    input:
        bam="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        bed="results/qc/rseqc/annotation.bed",
    threads: 32
    output:
        "results/qc/rseqc/{sample}_{unit}.readdistribution.txt",
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}_{unit}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log};
        """


rule rseqc_readdup:
    input:
        "results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readdup.DupRate_plot.pdf",
    threads: 32
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1;
        """

rule rseqc_readgc:
    input:
        "results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
    output:
        "results/qc/rseqc/{sample}_{unit}.readgc.GC_plot.pdf",
    threads: 32
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}_{unit}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_GC.py -i {input} -o {params.prefix} > {log} 2>&1;
        """

rule multiqc:
    input:
        get_multiqc_inputs,
    threads: 32
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.5.3/bio/multiqc"