rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}_{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}_{unit}.log",
    benchmark:
        "logs/star/{sample}_{unit}.bench.tsv",
    threads: 92
    params:
        idx=lambda wc, input: input.index,
        extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]} --runThreadN 92 ',
    wrapper:
        "v3.5.3/bio/star/align"
