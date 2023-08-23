# **ASPEN**

### Table of Contents
- [ASPEN - **A**tac **S**eq **P**ip**E**li**N**e](#aspen)
    - [1. Outline](#1-outline)
    - [2. Runtime details](#2-runtime-details)
      - [2.1 Load Module On Biowulf](#21-load-module-on-biowulf)
      - [2.2 Create Sample Manifest](#22-create-sample-manifest)
      - [2.3 Run Pipeline](#23-run-pipeline)

### 1. Outline

ASPEN or **A**tac **S**eq **P**ip**E**li**N**e is CCBR's pipeline to calls peaks for ATAC-Seq datasets. It currently accepts paired-end Illumina data and calls peak using MACS2 and Genrich peak callers. Below is a brief outline of the steps performed by the pipeline:

* Trim PE reads with CutAdapt
* Remove reads aligning to known blacklisted regions, if provided
* Align reads provided genome using bowtie2. This step generates multiple output files:
  * `tagAlign.gz`, which is a BED6 format file mainly required for MACS2 peak calling
  * `dedup.bam`, deduplicated BAM format file which may be required for downstream processing (eg. TOBIAS)
  * `qsorted.bam`, query sorted BAM file for Genrich peak calling
* Pre-peakcalling QC metrics:
  * FastQC is run pre- and post-trimming
  * Fragment length distribution is calculated using custom scripts
  * Preseq is run to estimate library complexity
* Post-peakcalling QC metrics:
  * TSS distributions are calculated for each replicate
  * FRiP is calculated for each replicate
  * FRiPextra calculations are performed if _fripextra_ config files are supplied
    * Fraction of reads in DHS regions
    * Fraction of reads in promoter regions
    * Fraction of reads in enhancer regions
* Peak calling: Peaks (NarrowPeak format) are called using MACS2 and Genrich. If multiple replicates exist per sample, consensus peaks are called (BED format).
* Peak annotation: ChIPseeker is used to annotate peaks if the genome is hg38/hg19/mm10
* Motif enrichment: Motif Enrichment is calculated using HOMER and AME (MEME suite)
* Report: MultiQC is used to generate a customized final HTML report

### 2. Runtime details

#### 2.1 Load module on BIOWULF

To clone the repo run:

```bash
% module load ccbrpipeliner
```

This will add `aspen` to your PATH environmental variable.

### 2.2 Create Sample Manifest

Once the data is stored on biowulf, sample manifest TSV (`samples.tsv`) can be created to have the following columns:

1. replicateName

2. sampleName: Multiple replicates can have the same sampleName

3. path_to_R1_fastq: absolute path is preferred

4. path_to_R2_fastq: **Required!** Pipeline currently **only** supports paired-end data

Note that:

* symlinks are created for R1 and R2 files from the sample manifest in the results folder. These symlinks have the filenames \<replicateName\>.R1.fastq.gz and \<replicateName\>.R2.fastq.gz, respectively. Thus, original filenames do not matter and original files do not need to be renamed.
* **replicateName** is used as prefix for individual peak calls
* **sampleName** is used as prefix for consensus peak calls

### 2.3 Run Pipeline

To get more information about how to run the pipeline simply `cd` to the above `CCBR_ATACseq` folder and run

```bash
% aspen --help

##########################################################################################

Welcome to
____ ____ ___  ____ _  _
|__| [__  |__] |___ |\ |
|  | ___] |    |___ | \|    v0.6.0

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
  bash ./aspen -w/--workdir=<WORKDIR> -m/--runmode=<RUNMODE>

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

--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder                    (--runmode=init only)
--useenvmod|-e  : use  option while running Snakemake. This is for using modules on HPC instead of containers(default).
--help|-h       : print this help

Example commands:
  bash ./aspen -w=/my/ouput/folder -m=init
  bash ./aspen -w=/my/ouput/folder -m=dryrun
  bash ./aspen -w=/my/ouput/folder -m=run

##########################################################################################

VersionInfo:
  python          : python/3.7
  snakemake       : snakemake
  pipeline_home   : /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/ASPEN/dev
  git commit/tag  : 2a86949a63101eae06d9f8b7e53aa0234bf83b00	v0.6.1-13-g2a86949
  aspen_version   : v0.6.0

##########################################################################################
```

`aspen` is a **Biowulf** and **FRCE** specific wrapper script to the pipeline. Essentially, to run the pipeline the user has to follow 3 steps:

1. **Initialize the output folder**:

   This can be done using the following command:

   ```bash
   % aspen -m=init -w=<path_to_output_folder>
   ```

   The above command will create `config.yaml` and `samples.tsv` in the output folder. Please edit these as per your requirements. You can replace the `sample.tsv` file in the output folder with the sample manifest created in the previous step outlined above.

2. **Dryrun**:

   To dry-run the pipeline, you can run the following command after initializing the output folder:

   ```bash
   % aspen -m=dryrun -w=<path_to_output_folder>
   ```

   This should list out the chain of jobs (DAG) that will be submitted to the job scheduler.

3. **RUN!!**:

   If the dry-run looks as expected, then you can submit the job using:

   ```bash
   % aspen -m=run -w=<path_to_output_folder>
   ```

This  will submit one  _master_ job to slurm, which will in turn  keep managing the entire pipeline and submit/monitor jobs to the job scheduler as and when required.

This snakemake pipeline is built to run on [Biowulf](https://hpc.nih.gov/) and [FRCE](https://ncifrederick.cancer.gov/staff/frce/). But, as it uses containers for all intermediate steps, the pipeline can be executed on any HPC with minimal edits to the config file.

For comments/suggestions/advice please reach out to [Vishal Koparde](mailto:vishal.koparde@nih.gov).

<hr>
<p align="center">
	<a href="#aspen">Back to Top</a>
</p>
