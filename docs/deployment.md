# Running ASPEN

To effectively run the ASPEN (ATAC-Seq PipEliNe) on the Biowulf High-Performance Computing (HPC) system, please follow the detailed user guide below:


## Prerequisites

- **Biowulf Account:** Ensure you have an active Biowulf account.
- **Data Preparation:** Store your raw ATAC-Seq paired-end FASTQ files in a directory accessible from Biowulf.

## Setting Up the Environment

### Load the ASPEN Module on Biowulf

To access ASPEN, load the `ccbrpipeliner` module:

```bash
module load ccbrpipeliner/7
```

This command adds aspen to your system's PATH, allowing you to execute pipeline commands directly.

> **Note**: If you're operating outside of Biowulf, ensure that dependencies such as snakemake, python, and singularity are installed and accessible in your system's PATH.

### Create a Sample Manifest
ASPEN requires a sample manifest file (`samples.tsv`) to identify and organize your input data. This tab-separated file should include the following columns:

- `replicateName`: Unique identifier for each replicate.
- `sampleName`: Identifier for the sample; multiple replicates can share the same sample name.
- `path_to_R1_fastq`: Absolute path to the Read 1 FASTQ file.
- `path_to_R2_fastq`: Absolute path to the Read 2 FASTQ file (required for paired-end data).

> **Note**: Symlinks for R1 and R2 files will be created in the results directory, named as <replicateName>.R1.fastq.gz and <replicateName>.R2.fastq.gz, respectively. Therefore, original filenames do not need to be altered.

> **Note**: The `replicateName` is used as a prefix for individual peak calls, while the `sampleName` serves as a prefix for consensus peak calls.

> **Note**: For differential ATAC analysis, create a `contrasts.tsv` file with two columns (Group1 and Group2 ... aka Sample1 and Sample2, without headers) and place it in the output directory after initialization. Ensure each group/sample in the contrast has at least two replicates, as DESeq2 requires this for accurate contrast calculations.

## Running the ASPEN Pipeline
ASPEN operates through a series of modes to facilitate various stages of the analysis.

### Initialize the Working Directory
Begin by initializing your working directory, which will house configuration files and results. Replace <path_to_output_folder> with your desired output directory path:

```bash
aspen -m=init -w=<path_to_output_folder>
```

This command generates a config.yaml and a placeholder `samples.tsv` in the specified directory. Edit these files to reflect your experimental setup, replacing the placeholder `samples.tsv` with your prepared manifest. If performing differential analysis, include the `contrasts.tsv` file at this stage.

> **Note**: To explore all possible options of the `aspen` command you can either run it without any arguments or run `aspen --help`

Here is what help looks like:
```bash 
##########################################################################################

Welcome to
____ ____ ___  ____ _  _
|__| [__  |__] |___ |\ |
|  | ___] |    |___ | \|    v1.0.6

A_TAC_S_eq A_nalysis P_ip_E_li_N_e

##########################################################################################

This pipeline was built by CCBR (https://bioinformatics.ccr.cancer.gov/ccbr)
Please contact Vishal Koparde for comments/questions (vishal.koparde@nih.gov)

##########################################################################################

Here is a list of genome supported by aspen:

  * hg19          [Human]
  * hg38          [Human]
  * mm10          [Mouse]
  * mmul10        [Macaca mulatta(Rhesus monkey) or rheMac10]
  * bosTau9       [Bos taurus(cattle)]

aspen calls peaks using the following tools:

 * MACS2
 * Genrich        [RECOMMENDED FOR USE]

USAGE:
  bash /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/ASPEN/.v1.0.6/aspen -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

Required Arguments:
1.  WORKDIR     : [Type: String]: Absolute or relative path to the output folder with write permissions.

2.  RUNMODE     : [Type: String] Valid options:
    * init      : initialize workdir
    * dryrun    : dry run snakemake to generate DAG
    * run       : run with slurm
    * runlocal  : run without submitting to sbatch
    ADVANCED RUNMODES (use with caution!!)
    * unlock    : unlock WORKDIR if locked by snakemake NEVER UNLOCK WORKDIR WHERE PIPELINE IS CURRENTLY RUNNING!
    * reconfig  : recreate config file in WORKDIR (debugging option) EDITS TO config.yaml WILL BE LOST!
    * reset     : DELETE workdir dir and re-init it (debugging option) EDITS TO ALL FILES IN WORKDIR WILL BE LOST!
    * printbinds: print singularity binds (paths)
    * local     : same as runlocal

Optional Arguments:

--genome|-g     : genome eg. hg38
--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder                    (--runmode=init only)
--help|-h       : print this help

Example commands:
  bash /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/ASPEN/.v1.0.6/aspen -w=/my/output/folder -m=init
  bash /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/ASPEN/.v1.0.6/aspen -w=/my/output/folder -m=dryrun
  bash /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/ASPEN/.v1.0.6/aspen -w=/my/output/folder -m=run

##########################################################################################

VersionInfo:
  python          : python/3.10
  snakemake       : snakemake
  pipeline_home   : /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/ASPEN/.v1.0.6
  git commit/tag  : f4366158ad972bc667422dcd4783bd69fa041556	v1.0.6
  aspen_version   : v1.0.6

##########################################################################################

```

### Dry Run the Pipeline
Before executing the full analysis, perform a dry run to visualize the workflow and identify potential issues:

```bash
aspen -m=dryrun -w=<path_to_output_folder>
```
This step outlines the sequence of tasks (Directed Acyclic Graph - DAG) without actual execution, allowing you to verify the planned operations.

### Execute the Pipeline
If the dry run output is satisfactory, proceed to run the pipeline:

```bash
aspen -m=run -w=<path_to_output_folder>
```
This command submits a master job to the Slurm workload manager, which orchestrates the entire analysis workflow, managing job submissions and monitoring progress.

- Optional Argument:

`--singcache` or `-c`: Specify a Singularity cache directory. The default is `/data/${USER}/.singularity` if available; otherwise, it defaults to `${WORKDIR}/snakemake/.singularity`.

**Example Command**:

```bash
aspen -m=run -w=/my/output/folder -c /data/${USER}/.singularity
```
This command runs the pipeline with the specified working directory and Singularity cache directory.

## Monitor ASPEN runs

To monitor the status of your ASPEN pipeline and its associated jobs on a Slurm-managed system, you can utilize the squeue and scontrol commands. The squeue command provides information about jobs in the scheduling queue, while scontrol offers detailed insights into specific jobs.

To view all your active and pending jobs, execute:

```bash
squeue -u $USER
```
This command lists all jobs submitted by your user account, displaying details such as job IDs, partitions, job names, user names, job states, and the nodes allocated.

For more granular information about a specific job, including its child jobs spawned by ASPEN, use the scontrol command:

```bash
scontrol show job <jobid>
```
Replace <jobid> with the specific Job ID of interest. This will provide comprehensive details about the job's configuration and status, aiding in effective monitoring and management of your ASPEN pipeline processes.

