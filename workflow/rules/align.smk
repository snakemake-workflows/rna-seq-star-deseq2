rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
    output:
        "results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "results/star/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        index=lambda wc, input: input.index,
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
            "resources/genome.gtf", config["params"]["star"]
        ),
    threads: 24
    wrapper:
        "0.77.0/bio/star/align"
