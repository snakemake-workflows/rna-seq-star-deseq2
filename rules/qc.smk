## RSEQC
rule rseqc_gtf2bed:
    input:
        config["ref"]["annotation"]
    output:
        "star/rseqc/annotation.bed"
    log:
        "logs/gtf2bed.log"
    conda:
        "../envs/bedops.yaml"
    shell:
        "gtf2bed < {input} > {output} 2> {log}"


rule rseqc_junction_annotation:
    input:
        bam="star/{sample}-{unit}/Aligned.out.bam",
        bed="star/rseqc/annotation.bed"
    output:
        'star/rseqc/{sample}-{unit}.junctionanno.junction.bed'
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.log"
    params:
        extra = r'-q 255',  # STAR uses 255 as a scrore for uniq mappers
        prefix= 'star/rseqc/{sample}-{unit}.junctionanno'
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log[0]} 2>&1"

        
rule rseqc_junction_saturation:
    input:
        bam="star/{sample}-{unit}/Aligned.out.bam",
        bed="star/rseqc/annotation.bed"
    output:
        'star/rseqc/{sample}-{unit}.junctionsat.junctionSaturation_plot.pdf'
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.log"
    params:
        extra = r'-q 255', 
        prefix = 'star/rseqc/{sample}-{unit}.junctionsat'
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} > {log} 2>&1"


rule rseqc_stat:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
    output:
        "star/rseqc/{sample}-{unit}.stats.txt"
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

        
rule rseqc_infer:
    input:
        bam="star/{sample}-{unit}/Aligned.out.bam",
        bed="star/rseqc/annotation.bed"
    output:
        "star/rseqc/{sample}-{unit}.infer_experiment.txt"
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"

        
rule rseqc_innerdis:
    input:
        bam="star/{sample}-{unit}/Aligned.out.bam",
        bed="star/rseqc/annotation.bed"
    output:
        "star/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt"
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}-{unit}.log"
    params:
        prefix = "star/rseqc/{sample}-{unit}.inner_distance_freq"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam="star/{sample}-{unit}/Aligned.out.bam",
        bed="star/rseqc/annotation.bed"
    output:
        "star/rseqc/{sample}-{unit}.readdistribution.txt"
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}-{unit}.log"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        "star/{sample}-{unit}/Aligned.out.bam"
    output:
        'star/rseqc/{sample}-{unit}.readdup.DupRate_plot.pdf'
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}-{unit}.log"
    params:
        prefix = "star/rseqc/{sample}-{unit}.readdup"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"

        
rule rseqc_readgc:
    input:
        "star/{sample}-{unit}/Aligned.out.bam"
    output:
        'star/rseqc/{sample}-{unit}.readgc.GC_plot.pdf'
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}-{unit}.log"
    params:
        prefix = "star/rseqc/{sample}-{unit}.readgc"
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"
        

rule multiqc:
    input:
        expand("star/{unit.sample}-{unit.unit}/Aligned.out.bam", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.junctionanno.junction.bed", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.junctionsat.junctionSaturation_plot.pdf", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.stats.txt", unit=units.itertuples()),
        expand("star/rseqc/{sample}-{unit}.inner_distance_freq.inner_distance.txt", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.readdup.DupRate_plot.pdf", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.readgc.GC_plot.pdf", unit=units.itertuples()),
        "logs/rseqc/rseqc_junction_annotation",
        "logs/rseqc/rseqc_innerdis"
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.log"
    shell
        "0.31.1/bio/multiqc"
