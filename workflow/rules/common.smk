import glob
import subprocess

import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)


def get_final_output():
    final_output = expand(
        "results/diffexp/{contrast}.diffexp.symbol.tsv",
        contrast=config["diffexp"]["contrasts"],
    )
    final_output.append("results/deseq2/normcounts.symbol.tsv")
    final_output.append("results/counts/all.symbol.tsv")
    final_output.append("results/qc/multiqc_report.html")

    if config["pca"]["activate"]:
        # get all the variables to plot a PCA for
        pca_variables = list(config["diffexp"]["variables_of_interest"])
        if config["diffexp"]["batch_effects"]:
            pca_variables.extend(config["diffexp"]["batch_effects"])
        if config["pca"]["labels"]:
            pca_variables.extend(config["pca"]["labels"])
        final_output.extend(
            expand("results/pca.{variable}.svg", variable=pca_variables)
        )
    if config["pathway"]["activate"]:
        final_output.append("results/pathway/hallmark_GSVA.csv")
        final_output.append("results/pathway/hallmark_GSVA.svg")
    return final_output


validate(samples, schema="../schemas/samples.schema.yaml")

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


wildcard_constraints:
    sample="|".join(samples["sample_name"]),
    unit="|".join(units["unit_name"]),


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample
        accession = unit["sra"]
        # Check if pe or se
        if is_sra_paired_end(wildcards.sample):
            # paired-end sample
            return expand(
                "sra/{accession}_{read}.fastq", accession=accession, read=[1, 2]
            )
        else:
            # single-end sample
            return "sra/{}.fastq".format(accession)

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0
    return files


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["fq2"].isnull()
    sra_null = sample_units["sra"].isnull()
    paired = ~fq2_null | ~sra_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


def is_sra_paired_end(sample):
    sample_units = units.loc[sample]
    
    sra_accession = sample_units["sra"]
    
    if sra_accession.isnull().all():
        raise ValueError(f"No SRA accession found for sample {sample}")
    
    # Dictionary to cache SRA accession and their LibraryLayout
    sra_layout_cache = {}

    def get_library_layout(sra):
        if sra in sra_layout_cache:
            return sra_layout_cache[sra]
        
        try:
            cmd = f'esearch -db sra -query "{sra}" | efetch -format runinfo | tail -n +2 | cut -d "," -f 16'
            library_layout = subprocess.check_output(cmd, shell=True).decode("utf-8").strip()
            
            if library_layout not in ["SINGLE", "PAIRED"]:
                raise ValueError(f"Unexpected LibraryLayout: {library_layout} for SRA accession {sra}")
            
            sra_layout_cache[sra] = library_layout
            return library_layout
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Failed to fetch LibraryLayout for {sra}: {e}")
    
    # Fetch the LibraryLayout for the given sample
    library_layout = get_library_layout(sra_accession.iloc[0])
    
    # Check for consistency among all SRA accessions for this sample
    for sra in sra_accession.dropna():
        if get_library_layout(sra) != library_layout:
            raise ValueError(f"Inconsistent LibraryLayout found among samples with SRA accessions.")
    
    return library_layout == "PAIRED"


def get_fq(wildcards):
    if config["trimming"]["activate"]:
        # activated trimming, use trimmed data
        if is_paired_end(wildcards.sample):
            # paired-end sample
            return dict(
                zip(
                    ["fq1", "fq2"],
                    expand(
                        "results/trimmed/{sample}_{unit}_{group}.fastq.gz",
                        group=["R1", "R2"],
                        **wildcards,
                    ),
                )
            )
        # single end sample
        return {
            "fq1": "results/trimmed/{sample}_{unit}_single.fastq.gz".format(**wildcards)
        }
    else:
        u = units.loc[(wildcards.sample, wildcards.unit)]
        if pd.isna(u["fq1"]):
            accession = u["sra"]
            if is_sra_paired_end(wildcards.sample):
                return dict(
                    zip(
                        ["fq1", "fq2"],
                        expand(
                            "sra/{accession}_{group}.fastq",
                            accession=accession,
                            group=["1", "2"],
                        ),
                    )
                )
            else:
                return {"fq1": f"sra/{accession}.fastq"}
        else:
            if not is_paired_end(wildcards.sample):
                return {"fq1": f"{u.fq1}"}
            else:
                return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}


def get_strandedness(units):
    if "strandedness" in units.columns:
        return units["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * units.shape[0]


def get_deseq2_threads(wildcards=None):
    # https://twitter.com/mikelove/status/918770188568363008
    few_coeffs = False if wildcards is None else len(get_contrast(wildcards)) < 10
    return 1 if len(samples) < 100 or few_coeffs else 6


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies


def get_fastqs(wc):
    if config["trimming"]["activate"]:
        return expand(
            "results/trimmed/{sample}/{unit}_{read}.fastq.gz",
            unit=units.loc[wc.sample, "unit_name"],
            sample=wc.sample,
            read=wc.read,
        )
    unit = units.loc[wc.sample]
    if all(pd.isna(unit["fq1"])):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand(
            "sra/{accession}_{read}.fastq", accession=accession, read=wc.read[-1]
        )
    fq = "fq{}".format(wc.read[-1])
    return units.loc[wc.sample, fq].tolist()


def get_contrast(wildcards):
    return config["diffexp"]["contrasts"][wildcards.contrast]
