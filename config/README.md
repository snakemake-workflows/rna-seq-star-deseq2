# General settings
To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample and unit sheet

* Add samples to `config/samples.tsv`. For each sample, the columns `sample_name`, and `condition` have to be defined. The `condition` (healthy/tumor, before Treatment / after Treatment) will be used as contrast for the DEG analysis in DESeq2. To include other relevant variables such as batches, add a new column to the sheet.
* For each sample, add one or more sequencing units (runs, lanes or replicates) to the unit sheet `config/units.tsv`. By activating or deactivating `mergeReads` in the `config/config.yaml`, you can decide wether to merge replicates or run them individually. For each unit, define adapters, and either one (column `fq1`) or two (columns `fq1`, `fq2`) FASTQ files (these can point to anywhere in your system). Alternatively, you can define an SRA (sequence read archive) accession (starting with e.g. ERR or SRR) by using a column `sra`. In the latter case, the pipeline will automatically download the corresponding paired end reads from SRA. If both local files and SRA accession are available, the local files will be preferred.
To choose the correct geneCounts produced by STAR, you can define the strandedness of a unit. STAR produces counts for unstranded ('None' - default), forward oriented ('yes') and reverse oriented ('reverse') protocols.  

Missing values can be specified by empty columns or by writing `NA`.

# DESeq scenario

To initialize the DEG analysis, you need to define a model in the `config/config.yaml`. The model can include all variables introduced as columns in `config/samples.tsv`.
* The standard model is `~condition` - to include a batch variable, write `~batch + condition`.

