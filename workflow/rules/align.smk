rule align:
    input:
        fq1=get_map_reads_input_R1,
        fq2=get_map_reads_input_R2,
        index="resources/star_genome"
    output:
        # see STAR manual for additional output files
        "results/star/{sample}/Aligned.out.bam",
        "results/star/{sample}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}.log"
    params:
        # path to STAR reference genome index
        index="resources/star_genome",
        # optional parameters
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
              "resources/genome.gtf", config["params"]["star"])
    threads: 12
    wrapper:
        "0.64.0/bio/star/align"
