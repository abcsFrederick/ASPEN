# ASPEN

**A**tac **S**eq **P**ip**E**li**N**e

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13755867.svg)](https://doi.org/10.5281/zenodo.13755867)
[![release](https://img.shields.io/github/v/release/CCBR/ASPEN?color=blue&label=latest%20release)](https://github.com/CCBR/ASPEN/releases/latest)

### Table of Contents

- [ASPEN](#aspen)
  - [1. Outline](#1-outline)
  - [2. Runtime details](#2-runtime-details)
    - [2.1 Load Module On Biowulf](#21-load-module-on-biowulf)
    - [2.2 Create Sample Manifest](#22-create-sample-manifest)
    - [2.3 Run Pipeline](#23-run-pipeline)
  - [3. Genomes](#3-genomes)
  - [4. Disclaimer](#4-disclaimer)
  - [5. Help](#5-help)

### 1. Outline

ASPEN or **A**tac **S**eq **P**ip**E**li**N**e is CCBR's pipeline to calls peaks for ATAC-Seq datasets. It currently accepts paired-end Illumina data and calls peak using [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137) and [Genrich](https://github.com/jsh58/Genrich) peak callers. Below is a brief outline of the steps performed by the pipeline:

- Trim PE reads with [CutAdapt](https://doi.org/10.14806/ej.17.1.200)
- Remove reads aligning to known [blacklisted regions](https://doi.org/10.1038/s41598-019-45839-z), if provided
- Align reads provided genome using [bowtie2](https://doi.org/10.1038%2Fnmeth.1923). This step generates multiple output files:
  - `tagAlign.gz`, which is a BED6 format file mainly required for MACS2 peak calling
  - `dedup.bam`, deduplicated BAM format file which may be required for downstream processing (eg. [TOBIAS](https://github.com/CCBR/CCBR_tobias))
  - `qsorted.bam`, query sorted BAM file for Genrich peak calling
- Pre-peakcalling QC metrics:
  - [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is run pre- and post-trimming
  - Fragment length distribution is calculated using custom scripts
  - [Preseq](https://smithlabresearch.org/software/preseq/) is run to estimate library complexity
- Post-peakcalling QC metrics:
  - TSS distributions are calculated for each replicate
  - FRiP (Fraction of Reads in Peaks) is calculated for each replicate
  - FRiPextra calculations are performed if _fripextra_ config files are supplied
    - Fraction of reads in DHS regions
    - Fraction of reads in promoter regions
    - Fraction of reads in enhancer regions
- Peak calling: Peaks (NarrowPeak format) are called using MACS2 and Genrich. If multiple replicates exist per sample, consensus peaks are called (BED format).
- Peak annotation: [ChIPseeker](https://doi.org/10.1093/bioinformatics/btv145) is used to annotate peaks if the genome is hg38/hg19/mm10
- Motif enrichment: Motif Enrichment is calculated using [HOMER](http://homer.ucsd.edu/homer/) and [AME (MEME suite)](https://meme-suite.org/meme/doc/ame.html)
- Report: [MultiQC](10.1093/bioinformatics/btw354) is used to generate a customized final HTML report

### 2. Runtime details

#### 2.1 Load module on BIOWULF

To clone the repo run:

```bash
% module load ccbrpipeliner
```

This will add `aspen` to your PATH environmental variable.

> NOTE:
> If not running on BIOWULF, please ensure that the required tools like snakemake, python and singularity are in PATH.

### 2.2 Create Sample Manifest

Once the data is stored on biowulf, sample manifest TSV (`samples.tsv`) can be created to have the following columns:

1. replicateName

2. sampleName: Multiple replicates can have the same sampleName

3. path_to_R1_fastq: absolute path is preferred

4. path_to_R2_fastq: **Required!** Pipeline currently **only** supports paired-end data

Note that:

- symlinks are created for R1 and R2 files from the sample manifest in the results folder. These symlinks have the filenames \<replicateName\>.R1.fastq.gz and \<replicateName\>.R2.fastq.gz, respectively. Thus, original filenames do not matter and original files do not need to be renamed.
- **replicateName** is used as prefix for individual peak calls
- **sampleName** is used as prefix for consensus peak calls

> NOTE:
> Optionally, if running differential ATAC please also provide `contrasts.tsv` in the output folder after running `init`. This is a simple tab-delimited text file with 2 columns (_Group1_ and _Group2_) without any headers.

### 2.3 Run Pipeline

To get more information about how to run the pipeline simply `cd` to the above `CCBR_ATACseq` folder and run

```bash
% aspen --help
##########################################################################################

Welcome to
____ ____ ___  ____ _  _
|__| [__  |__] |___ |\ |
|  | ___] |    |___ | \|    v1.0.0

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

--help|-h       : print this help
--genome|-g     : genome eg. hg38
--manifest|-s   : absolute path to samples.tsv. This will be copied to output folder                    (--runmode=init only)
--useenvmod|-e  : use "--use-enmodules" option while running Snakemake. This is for using modules on HPC instead of containers(default).
--singcache|-c  : singularity cache directory. Default is `/data/${USER}/.singularity` if available, or falls back to `${WORKDIR}/.singularity`.


Example commands:
  bash ./aspen -w=/my/output/folder -m=init
  bash ./aspen -w=/my/output/folder -m=dryrun
  bash ./aspen -w=/my/output/folder -m=run
  bash ./aspen -w=/my/output/folder -m=run -c /data/${USER}/.singularity

##########################################################################################

VersionInfo:
  python          : python/3.10
  snakemake       : snakemake
  pipeline_home   : /gpfs/gsfs10/users/CCBR_Pipeliner/Pipelines/ASPEN/dev
  git commit/tag  : 51cb3aee2142ce1226acca97e1e662d54c881a13    v1.0.0-2-g51cb3ae
  aspen_version   : v1.0.0

##########################################################################################
```

`aspen` is a **Biowulf** and **FRCE** specific wrapper script to the pipeline. Essentially, to run the pipeline the user has to follow 3 steps:

1. **Initialize the output folder**:

   This can be done using the following command:

   ```bash
   % aspen -m=init -w=<path_to_output_folder>
   ```

   The above command will create `config.yaml` and `samples.tsv` in the output folder. Please edit these as per your requirements. You can replace the `samples.tsv` file in the output folder with the sample manifest created in the previous step outlined above. `contrasts.tsv` should also be included if running differential ATAC.

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

This will submit one _master_ job to slurm, which will in turn keep managing the entire pipeline and submit/monitor jobs to the job scheduler as and when required.

### 3. Genomes

ASPEN supports the following genome versions:

| Genome Version | Organism                                  |
| -------------- | ----------------------------------------- |
| hg38           | Human                                     |
| hg19           | Human                                     |
| mm10           | Mouse                                     |
| mmul10         | Macaca mulatta(Rhesus monkey) or rheMac10 |
| bosTau9        | Bos taurus(cattle)                        |

### 4. Disclaimer

This snakemake pipeline is built to run on [Biowulf](https://hpc.nih.gov/) and [FRCE](https://ncifrederick.cancer.gov/staff/frce/). But, as it uses containers for all intermediate steps, the pipeline can be executed on any HPC with minimal edits to the config file.

### 5. Help

For comments/suggestions/advice please reach out to [Vishal Koparde](mailto:vishal.koparde@nih.gov) or [CCBR_Pipeliner](mailto:CCBR_Pipeliner@mail.nih.gov). You can also open a new issue [here](https://github.com/CCBR/ASPEN/issues).

<hr>
<p align="center">
	<a href="#aspen">Back to Top</a>
</p>

<hr>
