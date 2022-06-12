# **ASAP**

### Outline

ASAP or **A**tac **S**eq **A**nalysis **P**ipeline is CCBR's pipeline to calls peaks for ATACseq datasets. It currently accepts paired-end Illumina data and calls peak using MACS2 and Genrich peak callers. Below is a brief outline of the steps performed by the pipeline:

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

### Clone Repo

To clone the repo run:

```bash
% git clone git@github.com:kopardev/CCBR_ATACseq.git
```

This will create `CCBR_ATACseq` subfolder at your current location.

### Create Sample Manifest

Once the data is stored on biowulf, sample manifest TSV (`samples.tsv`) can be created to have the following columns:

1. replicateName

2. sampleName: Multiple replicates can have the same sampleName

3. path_to_R1_fastq: absolute path is preferred

4. path_to_R2_fastq: **Required!** Pipeline currently **only** supports paired-end data

Note that:

* symlinks are created for R1 and R2 files from the sample manifest in the results folder. These symlinks have the filenames \<replicateName\>.R1.fastq.gz and \<replicateName\>.R2.fastq.gz, respectively. Thus, original filenames do not matter and original files do not need to be renamed.
* **replicateName** is used as prefix for individual peak calls
* **sampleName** is used as prefix for consensus peak calls

### Run Pipeline

To get more information about how to run the pipeline simply `cd` to the above `CCBR_ATACseq` folder and run

```bash
% ./run_atac --help

Pipeline Dir: /gpfs/gsfs11/users/kopardevn/GitRepos/CCBR_ATACseq
Git Commit/Tag: 2ba69f785d0e2fd5c3f7c77de6807e0ef0d00312

#
#
#
  Unknown argument!
#
#
#

/gpfs/gsfs11/users/kopardevn/GitRepos/CCBR_ATACseq/run_atac
--> run CCBR ATACseq pipeline

USAGE:
  bash ./run_atac -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
Required Arguments:
1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) run : run with slurm
    *) reset : DELETE workdir dir and re-init it
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
    *) runlocal : run without submitting to sbatch
2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
```

`run_atac` is a **Biowulf** specific wrapper script to the pipeline. Essentially, to run the pipeline the user has to follow 3 steps:

1. **Initialize the output folder**:

   This can be done using the following command:

   ```bash
   % ./run_atac -m=init -w=<path_to_output_folder>
   ```

   The above command will create `config.yaml` and `samples.tsv` in the output folder. Please edit these as per your requirement. You can replace the `sample.tsv` file in the output folder with the sample manifest created in the previous step outlined above.

2. **Dryrun**:

   To dry-run the pipeline, you can run the following command after initializing the output folder:

   ```bash
   % ./run_atac -m=dryrun -w=<path_to_output_folder>
   ```

   This is should list out the chain of jobs that will be submitted to the job-scheduler.

3. **RUN!!**:

   If dry-run looks as expected, then you can submitted the job using:

   ```bash
   % ./run_atac -m=run -w=<path_to_output_folder>
   ```

This  will submit one  _master_ job to slurm, which will in-turn  keep managing the entire pipeline and submit/monitor jobs to the job scheduler as-and-when required.

### 
