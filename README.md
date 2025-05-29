# Snakemake workflow: rna-seq-star-deseq2 (using pcluster-slurm executor)

This is my adoption of the [original forked repo's rna seq start deseq2](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2) worflow, but to use AWS Parallel Cluster, via the [pcluser slurm snakemake executor](https://github.com/Daylily-Informatics/snakemake-executor-plugin-pcluster-slurm-ref).

This workflow performs a differential gene expression analysis with STAR and Deseq2.

# Usage
## Have an `AWS Parallel Cluster` ( using slurm as the scheduler ) Running.
From the headnode, proceed to the following steps.

## Prereq : Conda
- Have conda activated.

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

```bash
conda activate snakemake

# Set your cache dir for saving resources useful across other jobs, snakemake uses this when the `--cache` flag is set.

mkdir /fsx/resources/environments/containers/ubuntu/rnaseq_cache/
export SNAKEMAKE_OUTPUT_CACHE=/fsx/resources/environments/containers/ubuntu/rnaseq_cache/

# I set a partition relevant to my install, but if you specify nothing, you will get an error along the lines of <could not find appropriate nodes>.
snakemake --use-conda --use-singularity   --singularity-prefix /fsx/resources/environments/containers/ubuntu/ --singularity-args "  -B /tmp:/tmp -B /fsx:/fsx  -B /home/$USER:/home/$USER -B $PWD/:$PWD" --conda-prefix /fsx/resources/environments/containers/ubuntu/ --executor pcluster-slurm --default-resources slurm_partition=i64,i96,i192 --cache -p --verbose -k --max-threads 20000 --cores 20000 -j 14 -n   --conda-create-envs-only
```

- there seems to be a bug which requires you to run with  `--conda-create-envs-only` first...
- another bug with how snakemake detects max allowd threads per job limits the threads to the `nproc` of your head node.  Setting `--max-threads 20000 --cores 20000` gets around this crudely.


- Remove the `-n` flag, and run not in dryrun mode.
- `-j` sets the max jobs slurm will allow active at one time.
- Watch your running nodes/jobs using `squeue` (also, `q` cluster commands work, but not reliably and are not supported).
- **note:** it seems a bug in this example causes a few jobs to fail (investigating)

#### What Partitions Are Available?
Use `sinfo` to learn about your cluster (note, `sinfo` reports on all potential and active compute nodes. Read the docs to interpret which are active, which are not yet requested spot instances, etc). Below is what the [daylily AWS parallel cluster](https://github.com/Daylily-Informatics/daylily/blob/main/config/day_cluster/prod_cluster.yaml) looks like.

```bash
sinfo
PARTITION AVAIL  TIMELIMIT  NODES  STATE NODELIST
i8*          up   infinite     12  idle~ i8-dy-gb64-[1-12]
i64          up   infinite     16  idle~ i64-dy-gb256-[1-8],i64-dy-gb512-[1-8]
i96          up   infinite     16  idle~ i96-dy-gb384-[1-8],i96-dy-gb768-[1-8]
i128         up   infinite     28  idle~ i128-dy-gb256-[1-8],i128-dy-gb512-[1-10],i128-dy-gb1024-[1-10]
i192         up   infinite     30  idle~ i192-dy-gb384-[1-10],i192-dy-gb768-[1-10],i192-dy-gb1536-[1-10]
a192         up   infinite     30  idle~ a192-dy-gb384-[1-10],a192-dy-gb768-[1-10],a192-dy-gb1536-[1-10]
```
-  As I look at this, it is possible that if unset, the partition will default to `i8` in the output above. Maybe.

#### Budgets, and the `--comment` sbatch flag
I etensively make use of  [Cost allocation tags with AWS ParallelCluster](https://github.com/Daylily-Informatics/aws-parallelcluster-cost-allocation-tags) in the [daylily omics analysis framework](https://github.com/Daylily-Informatics/daylily?tab=readme-ov-file#daylily-aws-ephemeral-cluster-setup-0714) [_$3 30x WGS analysis_](https://github.com/Daylily-Informatics/daylily?tab=readme-ov-file#3-30x-fastq-bam-bamdeduplicated-snvvcfsvvcf-add-035-for-a-raft-of-qc-reports)  to track AWS cluster usage costs in realtime, and impose limits where appropriate (by user and project). This makes use of overriding the `--comment` flag to hold `project/budget` tags applied to ephemeral AWS resources, and thus enabling cost tracking/controls.

* To change the --comment flag in v`0.0.8` of the pcluster-slurm plugin, set the comment flag value in the envvar `SMK_SLURM_COMMENT=RandD` (RandD is the default).

 
