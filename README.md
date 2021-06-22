# **CCBR ATACseq Pipeline**

### Clone Repo

To clone the repo run:

```bash
% git clone git@github.com:kopardev/CCBR_ATACseq.git
```

This will create `CCBR_ATACseq` subfolder at your current location.

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

   The above command will create `config.yaml` and `samples.tsv` in the output folder. Please edit these as per your requirement. `samples.tsv` is automatically populated with sample dataset, which can be replaced with your data. The 4 columns in the file are:

   1. sampleName: This is name of the sample. In case of replicates, this is the replicate name.
   2. groupName: If no replicates are present, then this is same as sampleName, else this is actually the name of the sample represented by multiple replicates
   3. path_to_R1_fastq: absolute path is preferred
   4. path_to_R2_fastq: Required! Pipeline currently only supports paired-end data

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