rule get_genome:
    output:
        "resources/genome.fasta",
    log:
        "logs/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: "omit-software"
    wrapper:
        "v7.2.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: "omit-software"
    log:
        "logs/get_annotation.log",
    wrapper:
        "v7.2.0/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    wrapper:
        "v7.2.0/bio/samtools/faidx"


rule bwa_index:
    input:
        "resources/genome.fasta",
    output:
        multiext("resources/genome.fasta", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v7.2.0/bio/bwa/index"


rule star_index:
    input:
        fasta="resources/genome.fasta",
        gtf="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    log:
        "logs/star_index_genome.log",
    cache: True
    params:
        extra=lookup(within=config, dpath="params/star/index", default=""),
    threads: 4
    wrapper:
        "v7.2.0/bio/star/index"
