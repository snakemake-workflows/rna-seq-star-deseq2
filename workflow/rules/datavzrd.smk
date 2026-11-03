rule datavzrd_diffexp:
    input:
        config=workflow.source_path(
            "../resources/datavzrd/diffexp-template.datavzrd.yaml"
        ),
        table="results/diffexp/{contrast}.diffexp.symbol.tsv",
    output:
        report(
            directory("results/datavzrd-report/{contrast}.diffexp"),
            htmlindex="index.html",
            caption="../report/diffexp-datavzrd.rst",
            category="Differential expression",
            labels=lambda wildcards: {"contrast": wildcards.contrast},
        ),
    log:
        "logs/datavzrd/{contrast}.diffexp.log",
    wrapper:
        "master/utils/datavzrd"
