rule align:
    input:
        unpack(get_fq),
        idx="resources/star_genome",
    output:
        aln="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]}',
    threads: 24
    wrapper:
        "v7.2.0/bio/star/align"
