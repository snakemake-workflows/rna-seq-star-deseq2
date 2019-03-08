## RSEQC
rule rseqc_gtf2bed12:
    input:
        config["ref"]["annotation"]
    output:
        "star/rseqc/annotation.bed"
    log:
        "logs/gtf2bed12.stderr"
    params:
        script = "../scripts/gtf2bed.pl"
    shell:
        """
        perl {params.script} {input} > {output} 2> {log}
        """


rule rseqc_junction_annotation:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/rseqc/annotation.bed"
    output:
        touch('star/rseqc/{sample}-{unit}.junction_annotation.done')
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.stdout",
        "logs/rseqc/rseqc_junction_annotation/{sample}-{unit}.stderr"
    params:
        extra = r'-q 255',  # STAR uses 255 as a scrore for uniq mappers
        prefix= 'star/rseqc/{sample}-{unit}.junctionanno'
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        junction_annotation.py {params.extra} -i {input[0]} -r {input[1]} -o {params.prefix} > {log[0]} 2> {log[1]};
        """

        
rule rseqc_junction_saturation:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/rseqc/annotation.bed"
    output:
        touch('star/rseqc/{sample}-{unit}.junction_saturation.done')
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.stdout",
        "logs/rseqc/rseqc_junction_saturation/{sample}-{unit}.stderr"
    params:
        extra = r'-q 255', 
        prefix = 'star/rseqc/{sample}-{unit}.junctionsat'
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        junction_saturation.py {params.extra} -i {input[0]} -r {input[1]} -o {params.prefix} > {log[0]} 2> {log[1]};
        """


rule rseqc_stat:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
    output:
        "star/rseqc/{sample}-{unit}.stats.txt"
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}-{unit}.stderr"
    params:
        extra = r''
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        bam_stat.py {params.extra} -i {input} > {output} 2> {log}
        """

        
rule rseqc_infer:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/rseqc/annotation.bed"
    output:
        "star/rseqc/{sample}-{unit}.infer_experiment.txt"
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}-{unit}.stderr"
    params:
        extra = r''
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        infer_experiment.py {params.extra} -r {input[1]} -i {input[0]} > {output} 2> {log}
        """

        
rule rseqc_genebody:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/rseqc/annotation.bed"
    output:
        touch("star/rseqc/{sample}-{unit}.genebody.done")
    priority: 1
    log:
        "logs/rseqc/rseqc_genebody/{sample}-{unit}.stdout",
        "logs/rseqc/rseqc_genebody/{sample}-{unit}.stderr"
    params:
        extra = r'',
        prefix = "star/rseqc/{sample}-{unit}.geneBodyCoverage"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        geneBody_coverage.py {params.extra} -r {input[1]} -i {input[0]} -o {params.prefix} > {log[0]} 2> {log[1]}
        """


rule rseqc_innerdis:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/rseqc/annotation.bed"
    output:
        touch('star/rseqc/{sample}-{unit}.innerdistance.done')
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}-{unit}.stdout",
        "logs/rseqc/rseqc_innerdis/{sample}-{unit}.stderr"
    params:
        extra = r'',
        prefix = "star/rseqc/{sample}-{unit}.inner_distance_freq"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        inner_distance.py {params.extra} -r {input[1]} -i {input[0]} -o {params.prefix} > {log[0]} 2> {log[1]}
        """


rule rseqc_readdis:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/rseqc/annotation.bed"
    output:
        "star/rseqc/{sample}-{unit}.readdistribution.txt"
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}-{unit}.stderr"
    params:
        extra = r''
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_distribution.py {params.extra} -r {input[1]} -i {input[0]} > {output} 2> {log}
        """


rule rseqc_readdup:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
    output:
        touch('star/rseqc/{sample}-{unit}.readdup.done')
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}-{unit}.stdout",
        "logs/rseqc/rseqc_readdup/{sample}-{unit}.stderr"
    params:
        extra = r'',
        prefix = "star/rseqc/{sample}-{unit}.readdup"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_duplication.py {params.extra} -i {input[0]} -o {output} > {log[0]} 2> {log[1]}
        """

        
rule rseqc_readgc:
    input:
        "star/{sample}-{unit}/Aligned.out.bam",
    output:
        touch('star/rseqc/{sample}-{unit}.readgc.done')
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}-{unit}.stdout",
        "logs/rseqc/rseqc_readgc/{sample}-{unit}.stderr"
    params:
        extra = r'',
        prefix = "star/rseqc/{sample}-{unit}.readgc"
    conda:
        "../envs/rseqc.yaml"
    shell:
        """
        read_GC.py {params.extra} -i {input[0]} -o {params.prefix} > {log[0]} 2> {log[1]}
        """
        

rule multiqc:
    input:
        expand("star/{unit.sample}-{unit.unit}/Log.final.out", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.junction_annotation.done", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.junction_saturation.done", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.infer_experiment.txt", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.stats.txt", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.innerdistance.done", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.readdistribution.txt", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.readdup.done", unit=units.itertuples()),
        expand("star/rseqc/{unit.sample}-{unit.unit}.readgc.done", unit=units.itertuples())
    output:
        "qc/multiqc_report.html"
    log:
        "logs/multiqc.stdout",
        "logs/multiqc.stderr"
    conda:
        "../envs/multiqc.yaml"
    params:
        extra = '',
        out = "qc",
        align = "star",
        rseqc = "star/rseqc",
        rseqc_anno = "logs/rseqc/rseqc_junction_annotation",
        rseqc_dis = "logs/rseqc/rseqc_innerdis"
    shell:
        """
        multiqc {params.extra} -f -o {params.out} {params.align} {params.rseqc_anno} {params.rseqc_dis} {params.rseqc} > {log[0]} 2> {log[1]}
        """
