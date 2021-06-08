rule align_pe:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        index="resources/star_genome",
    output:
        "results/star/pe/{sample}-{unit}/Aligned.out.bam",
        "results/star/pe/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-pe/{sample}-{unit}.log",
    params:
        index=lambda wc, input: input.index,
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.64.0/bio/star/align"


rule align_se:
    input:
        fq1=get_map_reads_input_R1,
        index="resources/star_genome",
    output:
        "results/star/se/{sample}-{unit}/Aligned.out.bam",
        "results/star/se/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star-se/{sample}-{unit}.log",
    params:
        index=lambda wc, input: input.index,
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.64.0/bio/star/align"
