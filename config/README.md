# General configuration

To configure this workflow, modify `config/config.yaml` according to your needs, following the explanations provided in the file.

## `DESeq2` differential expression analysis configuration

To successfully run the differential expression analysis, you will need to tell DESeq2 which sample annotations to use (annotations are columns in the `samples.tsv` file described below).
This is done in the `config.yaml` file with the entries under `diffexp:`.
The comments for the entries should give all the necessary infos and linkouts.
But if in doubt, please also consult the [`DESeq2` manual](https://www.bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

# Sample and unit setup

The sample and unit setup is specified via tab-separated tabular files (`.tsv`).
Missing values can be specified by empty columns or by writing `NA`.

## sample sheet

The default sample sheet is `config/samples.tsv` (as configured in `config/config.yaml`).
Each sample refers to an actual physical sample, and replicates (both biological and technical) may be specified as separate samples.
For each sample, you will always have to specify a `sample_name`.
In addition, all `variables_of_interest` and `batch_effects` specified in the `config/config.yaml` under the `diffexp:` entry, will have to have corresponding columns in the `config/samples.tsv`.
Finally, the sample sheet can contain any number of additional columns.
So if in doubt about whether you might at some point need some metadata you already have at hand, just put it into the sample sheet already---your future self will thank you.

## unit sheet

The default unit sheet is `config/units.tsv` (as configured in `config/config.yaml`).
For each sample, add one or more sequencing units (for example if you have several runs or lanes per sample).

### `.fastq` file source

For each unit, you will have to define a source for your `.fastq` files.
This can be done via the columns `fq1`, `fq2` and `sra`, with either of:
1. A single `.fastq` file for single-end reads (`fq1` column only; `fq2` and `sra` columns present, but empty).
  The entry can be any path on your system, but we suggest something like a `raw/` data directory within your analysis directory.
2. Two `.fastq` files for paired-end reads (columns `fq1` and `fq2`; column `sra` present, but empty).
  As for the `fq1` column, the `fq2` column can also point to anywhere on your system.
3. A sequence read archive (SRA) accession number (`sra` column only; `fq1` and `fq2` columns present, but empty).
  The workflow will automatically download the corresponding `.fastq` data (currently assumed to be paired-end).
  The accession numbers usually start with SRR or ERR and you can find accession numbers for studies of interest with the [SRA Run Selector](https://trace.ncbi.nlm.nih.gov/Traces/study/).
If both local files and an SRA accession are specified for the same unit, the local files will be used.

### strandedness of library preparation protocol

To get the correct `geneCounts` from `STAR` output, you can provide information on the strandedness of the library preparation protocol used for a unit.
`STAR` can produce counts for unstranded (`none` - this is the default), forward oriented (`yes`) and reverse oriented (`reverse`) protocols.  
Enter the respective value into a `strandedness` column in the `units.tsv` file.

### adapter trimming and read filtering

Finally, you can provide settings for the adapter trimming with `fastp` (see the [`fastp` documentation](https://github.com/OpenGene/fastp)) via the `units.tsv` columns `fastp_adapters` and `fastp_extra`.
However, if you leave those two columns empty (no whitespace!), `fastp` will auto-detect adapters and the workflow will set sensible defaults for trimming of RNA-seq data.
If you use this automatic inference, make sure to double-check the `Detected read[12] adapter:` entries in the resulting `fastp` HTML report.
This is part of the final `snakemake` report of the workflow, or can be found in the sample-specific folders under `results/trimmed/`, once a sample has been processed this far.
If the auto-detection didn't work at all (empty `Detected read[12] adapter:` entries), or the `Occurrences` in the `Adapters` section are lower than you would expect, please ensure that you find out which adapters were used and configure the adapter trimming manually:

In the column `fastp_adapters`, you can specify [known adapter sequences to be trimmed off by `fastp`](https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters), including the command-line argument for the trimming.
For example, specify the following string in this column: `--adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`
If you don't know the adapters used, leave this empty (an empty string, containing no whitespace), and `fastp` will auto-detect the adapters that need to be trimmed.
If you want to make the auto-detection explicit for paired-end samples, you can also specify `--detect_adapter_for_pe`.

In the column `fastp_extra`, you can specify [further `fastp` command-line settings](https://github.com/OpenGene/fastp?tab=readme-ov-file#all-options).
If you leave this empty (an empty string, containing no whitespace), the workflow will set its default:
* [`--trim_poly_x --poly_x_min_len 7`](https://github.com/OpenGene/fastp?tab=readme-ov-file#polyx-tail-trimming): This poly-X trimming removes polyA tails if they are 7 nucleotides or longer.
  It is run after adapter trimming.
* [`--trim_poly_g --poly_g_min_len 7`](https://github.com/OpenGene/fastp?tab=readme-ov-file#polyx-tail-trimming): This poly-G trimming removes erroneous G basecalls at the tails of reads, if they are 7 nucleotides or longer.
  These Gs are artifacts in Illumina data from [machines with a one channel or two channel color chemistry](https://github.com/OpenGene/fastp/pull/508#issuecomment-3028078859).
  We currently set this by default, because [the auto-detection for the respective machines is lacking the latest machine types](https://github.com/OpenGene/fastp/pull/508).
  When the above-linked pull request is updated and merged, we can remove this and rely on the auto-detection.
If you want to specify additional command line options, we recommend always including those parameters in your units.tsv, as well.
Here's the full concatenation for copy-pasting:

```bash
--trim_poly_x --poly_x_min_len 7 --trim_poly_g --poly_g_min_len 7
```

#### Lexogen 3' QuantSeq adapter trimming

For this data, adapter trimming should automatically work as expected with the use of `fastp`.
The above-listed defaults are equivalent to an adaptation of the [Lexogen read preprocessing recommendations for 3' FWD QuantSeq data with `cutadapt`](https://faqs.lexogen.com/faq/what-sequences-should-be-trimmed).
The only difference is that we don't do any length filtering with these defaults.
If you want to exactly mirror the Lexogen recommendations, please use this for the `fastp_extra` column in your `units.tsv`:

```bash
--length_required 20 --trim_poly_x --poly_x_min_len 7 --trim_poly_g --poly_g_min_len 7
```

The `fastp` equivalents, including minimal deviations from the recommendations, are motivated as follows:
* `-m`: In cutadapt, this is the short version of `--minimum-length`. The `fastp` equivalent is `--length_required`.
* `-O`: Here, `fastp` doesn't have an equivalent option, so we currently have to live with the suboptimal default of `4`. This is greater than the `min_overlap=3` used here; but smaller than the value of `7`, a threshold that we have found avoids removing randomly matching sequences when combined with the typical Illumina `max_error_rate=0.005`.
* `-a "polyA=A{20}"`: This can be replaced by `fastp`'s dedicated option for `--trim_poly_x` tail removal ([which is run after adapter trimming](https://github.com/OpenGene/fastp?tab=readme-ov-file#global-trimming)).
* `-a "QUALITY=G{20}"`: This can be replaced by `fastp`'s dedicated option for the removal of artifactual trailing `G`s in Illumina data from [machines with a one channel or two channel color chemistry](https://github.com/OpenGene/fastp/pull/508#issuecomment-3028078859): `--trim_poly_g`.
  This is automatically activated for earlier Illumina machine models with this chemistry, but we recommend to activate it manually in the `fastp_extra` column of your `config/units.tsv` file for now, as [there are newer models that are not auto-detected, yet](https://github.com/OpenGene/fastp/pull/508).
  Also, we recommend to set `--poly_g_min_len 7`, to avoid trimming spurious matches of G-only stretches at the end of reads.
* `-n`: With the dedicated `fastp` options [getting applied in the right order](https://github.com/OpenGene/fastp?tab=readme-ov-file#global-trimming), this option is not needed any more.
* `-a "r1adapter=A{18}AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=3;max_error_rate=0.100000"`: We remove A{18}, as this is handled by `--trim_poly_x`.
  `fastp` uses the slightly higher `min_overlap` equivalent of `4`, which is currently hard-coded (and not exposed as a command-line argument).
  Because of this, we cannot set the `max_error_rate` to the Illumina error rate of about `0.005`.
* `-g "r1adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC;min_overlap=20"`: This is not needed any more, as `fastp` searches the read sequence for adapter sequences from the start of the read (see [the `fastp` adapter search code](https://github.com/OpenGene/fastp/blob/723a4293a42f1ca05b93c37db6a157b4235c4dcc/src/adaptertrimmer.cpp#L92)).
* `--discard-trimmed`: We omit this, as adapter sequence removal early in the read will leave short remaining read sequences that are subsequently filtered by `--length_required`.

