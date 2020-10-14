def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        u = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
        if is_single_end(**wildcards):
            return { 'fq1': f"{u.fq1}" }
        else:
            return { 'fq1': f"{u.fq1}",
                     'fq2': f"{u.fq2}" }

    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return dict(zip(
                ['fq1', 'fq2' ],
                expand("trimmed/{sample}-{unit}.{group}.fastq.gz", group=[1, 2], **wildcards)))
        # single end sample
        return { 'fq1': "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards) }


rule align:
    input:
        unpack(get_fq),
        index=config["ref"]["index"]
    output:
        # see STAR manual for additional output files
        "star/{sample}-{unit}/Aligned.sortedByCoord.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: 24
    wrapper:
        "0.66.0/bio/star/align"
