# Snakemake workflow: rna-seq-star-deseq2 (using pcluster-slurm executor) `0.1.2`

This is my adoption of the [original forked repo's rna seq start deseq2](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2) worflow, but to use AWS Parallel Cluster, via the [pcluser slurm snakemake executor](https://github.com/Daylily-Informatics/snakemake-executor-plugin-pcluster-slurm-ref).

This workflow performs a differential gene expression analysis with STAR and Deseq2.

# Prerequisites
## Have an `AWS Parallel Cluster` ( using slurm as the scheduler ) Running.

### _from scratch_ Use AWS Parallel Cluster (less reccomended)
- [AWS PC](https://aws.amazon.com/hpc/parallelcluster/)

### _pre-configured_ Use `daylily-ephemeral-cluster` (reccomended)

- Follow the setup instructions found here: [https://github.com/Daylily-Informatics/daylily-ephemeral-cluster](https://github.com/Daylily-Informatics/daylily-ephemeral-cluster).


## Conda
- If you have installed daylily-ehpemeral-cluster, once you log into the headnode, conda should be activated.
- If you roll your own, you'll need to install miniconda, and activate.

# Usage

## Clone Repo (it includes sample data)

```bash
git clone git@github.com:Daylily-Informatics/rna-seq-star-deseq2.git
cd rna-seq-star-deseq2
```

## Build The Snakemake (v8.*) Conda Env
```bash
conda create -n snakemake -c conda-forge  snakemake==8.24 snakedeploy tabulate yaml
conda activate snakemake
pip install snakemake-executor-plugin-pcluster-slurm==0.0.25

conda activate snakemake
snakemake --version
8.24  # was 8.20.6
```

### Run Test Data Workflow
_you are advised to run the following in a tmux or screen session_


#### Prepare Cache and TMPDIR

```bash
conda activate snakemake

# Set your cache dir for saving resources useful across other jobs, snakemake uses this when the `--cache` flag is set.

mkdir /fsx/resources/environments/containers/ubuntu/rnaseq_cache/
export SNAKEMAKE_OUTPUT_CACHE=/fsx/resources/environments/containers/ubuntu/rnaseq_cache/
export TMPDIR=/fsx/scratch/
```

#### Prepare `units.tsv`

```bash
cp config/units.tsv.template config/units.tsv
[[ "$(uname)" == "Darwin" ]] && sed -i "" "s|REGEX_PWD|$PWD|g" config/units.tsv || sed -i "s|REGEX_PWD|$PWD|g" config/units.tsv
```

#### Build Conda Env Caches 
_this can take ~1hr_

```bash
# I set partitions relevant to my AWS parallel cluster, but if you specify nothing, you will get an error along the lines of <could not find appropriate nodes>.
snakemake --use-conda --use-singularity   --singularity-prefix /fsx/resources/environments/containers/ubuntu/ --singularity-args "  -B /tmp:/tmp -B /fsx:/fsx  -B /home/$USER:/home/$USER -B $PWD/:$PWD" --conda-prefix /fsx/resources/environments/containers/ubuntu/ --executor pcluster-slurm --default-resources slurm_partition=i128,i192 --cache -p --verbose -k --max-threads 20000 --cores 20000 -j 14 -n   --conda-create-envs-only
```

- there seems to be a bug which requires you to run with  `--conda-create-envs-only` first, then once all envs are built, run the command.
- another bug with how snakemake detects max allowd threads per job limits the threads to the `nproc` of your head node.  Setting `--max-threads 20000 --cores 20000` gets around this crudely.

#### Prepare To Run The Command

- Remove the `-n` flag, and run not in dryrun mode.
- `-j` sets the max jobs slurm will allow active at one time.
- Watch your running nodes/jobs using `squeue` (also, `q` cluster commands work, but not reliably and are not supported).

##### What Partitions Are Available?
Use `sinfo` to learn about your cluster (note, `sinfo` reports on all potential and active compute nodes. Read the docs to interpret which are active, which are not yet requested spot instances, etc). Below is what the [daylily AWS parallel cluster](https://github.com/Daylily-Informatics/daylily/blob/main/config/day_cluster/prod_cluster.yaml) looks like.

```bash
sinfo
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
i8*          up   infinite     12  idle~ i8-dy-gb64-[1-12]
i64          up   infinite     16  idle~ i64-dy-gb256-[1-8],i64-dy-gb512-[1-8]
i128         up   infinite     28  idle~ i128-dy-gb256-[1-8],i128-dy-gb512-[1-10],i128-dy-gb1024-[1-10]
i192         up   infinite     30  idle~ i192-dy-gb384-[1-10],i192-dy-gb768-[1-10],i192-dy-gb1536-[1-10]
```

- Use the strings in `PARTITION`, ie: `i192` in the `slurm_partition=` config passed to snakemake.

##### Budgets, and the `--comment` sbatch flag
`daylily` makes extensive use of  [Cost allocation tags with AWS ParallelCluster](https://github.com/Daylily-Informatics/aws-parallelcluster-cost-allocation-tags) in the [daylily omics analysis framework](https://github.com/Daylily-Informatics/daylily?tab=readme-ov-file#daylily-aws-ephemeral-cluster-setup-0714) [_$3 30x WGS analysis_](https://github.com/Daylily-Informatics/daylily?tab=readme-ov-file#3-30x-fastq-bam-bamdeduplicated-snvvcfsvvcf-add-035-for-a-raft-of-qc-reports)  to track AWS cluster usage costs in realtime, and impose limits where appropriate (by user and project). This makes use of overriding the `--comment` flag to hold `project/budget` tags applied to ephemeral AWS resources, and thus enabling cost tracking/controls.

* To change the --comment flag in v`0.0.8` of the pcluster-slurm plugin, set the comment flag value in the envvar `SMK_SLURM_COMMENT=RandD` (RandD is the default).

##### Run The Command

```bash
snakemake --use-conda --use-singularity   --singularity-prefix /fsx/resources/environments/containers/ubuntu/ --singularity-args "  -B /tmp:/tmp -B /fsx:/fsx  -B /home/$USER:/home/$USER -B $PWD/:$PWD" --conda-prefix /fsx/resources/environments/containers/ubuntu/ --executor pcluster-slurm --default-resources slurm_partition=i128,i192 --cache -p --verbose -k --max-threads 20000 --cores 20000 -j 14 
```

 - You can watch progress with `watch squeue`.

#### Run w/Your Data

- Update the `config/units.tsv` (holds sample data location and other details) and `config/samples.tsv` (holds sample annotations).
- Edit `config/config.yaml` to change aspects of the pipeline.
- Run the `snakemake` command above *with* `-n`, tweak `-j` as needed, and if all looks good, run w/out `-n`.
