rule arriba:
    input:
        bam="results/star/{sample}/Aligned.out.bam",
        genome="resources/genome.fasta",
        annotation="resources/genome.gtf"
    output:
        fusions="results/arriba/{sample}.fusions.tsv",
        discarded="results/arriba/{sample}.fusions.discarded.tsv"
    params:
        blacklist=config["fusion"]["arriba"]["blacklist"],
        extra=config["params"]["fusion"]["arriba"]
    log:
        "results/logs/arriba/{sample}.log"
    threads: 1
    wrapper:
        "0.64.1/bio/arriba"
        
## TODO: Update
# rule fusioncatcher:
#    input:
#        fq1=lambda w: units.loc[(w.sample, "RNA"), "fq1"],
#        fq2=lambda w: units.loc[(w.sample, "RNA"), "fq2"]
#    output:
#        directory("fusioncatcher/{sample}")
#    params:
#        extra="-T tmp -d ../../fusioncatcher_data"
#    log:
#        "logs/fusioncatcher/{sample}.log"
#    threads:
#        8
#    shell:
#        "fusioncatcher -i {input.fq1},{input.fq2} -o {output} {params.extra} -p {threads} > {log}"