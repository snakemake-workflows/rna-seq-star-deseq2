rule align:
    input:
        unpack(get_fq),
        index="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        idx=lambda wc, input: input.index,
        extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]}',
    threads: 24
    wrapper:
        "v3.5.3/bio/star/align"
