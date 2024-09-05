rule GSVA:
    input:
        "results/deseq2/normcounts.tsv"
    output:
        "results/pathway/hallmark_GSVA.csv",
        "results/pathway/hallmark_GSVA.svg",
    log:
        "logs/pathway/gsva.log",
    conda:
        "../envs/pathway.yaml"
    script:
        "../scripts/pathway.R"