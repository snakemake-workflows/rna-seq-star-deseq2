rule star_align:
    input:
        unpack(get_fq),
        idx="resources/star_genome",
        gtf="resources/genome.gtf",
    output:
        aln="results/star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        reads_per_gene="results/star/{sample}-{unit}/ReadsPerGene.out.tab",
    log:
        "logs/star/{sample}-{unit}.log",
    params:
        extra=lambda wc, input: " ".join(
            [
                "--outSAMtype BAM SortedByCoordinate",
                "--quantMode GeneCounts",
                f'--sjdbGTFfile "{input.gtf}"',
                lookup(within=config, dpath="params/star/align", default=""),
            ]
        ),
    threads: 24
    wrapper:
        "v7.2.0/bio/star/align"
