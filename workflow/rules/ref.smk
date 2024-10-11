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
    cache: True
    threads: 32
    wrapper:
        "v3.5.3/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/genome.gtf",
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    threads: 32
    log:
        "logs/get_annotation.log",
    wrapper:
        "v3.5.3/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        "resources/genome.fasta",
    output:
        "resources/genome.fasta.fai",
    log:
        "logs/genome-faidx.log",
    cache: True
    threads: 32
    wrapper:
        "v3.5.3/bio/samtools/faidx"


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
    threads: 32
    wrapper:
        "v3.5.3/bio/bwa/index"


rule star_index:
    input:
        fasta="resources/genome.fasta",
        annotation="resources/genome.gtf",
    output:
        directory("resources/star_genome"),
    params:
        extra=lambda wc, input: f"--sjdbGTFfile {input.annotation} --sjdbOverhang 100",
    log:
        "logs/star_index_genome.log",
    cache: True
    threads: 32
    wrapper:
        "v3.5.3/bio/star/index"
