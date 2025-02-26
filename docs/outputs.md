# ASPEN Outputs

## Workdir

The workdir which is supplied as `-w` while running aspen `init`, `dryrun` and `run` commands will contain the following files:

```bash
WORKDIR
├── cluster.json
├── config.yaml
├── contrasts.tsv
├── dryrun_git_commit.txt
├── dryrun.log
├── fastqs
├── logs
├── results
├── runinfo.yaml
├── sampleinfo.txt
├── samples.tsv
├── scripts
├── slurm-48768339.out
├── snakemake.log
├── snakemake.stats
├── stats
├── submit_script.sbatch
└── tools.yaml
```

Here are more details about these files:

| File | File Type | Mode (`-m`) When This File is Created/Overwritten | Description |
|------|----------|--------------------------------------|-------------|
| cluster.json | JSON | init | Defines cluster resources per snakemake rule |
| config.yaml | YAML | init; can be edited later | Configurable parameters for this specific run |
| contrasts.tsv | TSV | Needs to be added in after init | List of contrasts to run, one per line; has no header |
| dryrun_git_commit.txt | TXT | dryrun | The git commit hash of the version of ASPEN used at dryrun |
| dryrun.log | TXT | dryrun | Log from `-m=dryrun` |
| fastqs | FOLDER | dryrun | Folder containing symlinks to raw data |
| logs | FOLDER | dryrun | Folder containing all logs including Slurm `.out` and `.err` files |
| results | FOLDER | Created at dryrun but populated during run | Main outputs folder |
| runinfo.yaml | YAML | After completion of run | Metadata about the run executor, etc. |
| sampleinfo.txt | TXT | dryrun, run | Tab-delimited mappings between `replicateNames` and `sampleNames` |
| samples.tsv | TSV | init; can be edited later | Tab-delimited manifest with `replicateName`, `sampleName`, `path_to_R1_fastq`, `path_to_R2_fastq`. This file has a header. |
| scripts | FOLDER | init | Folder keeps local copy of scripts called by various rules |
| slurm-49051815.out | TXT | run | Slurm `.out` file for the master job |
| snakemake.log | TXT | run | Snakemake `.log` file for the master job; older copies timestamped and moved into `logs` folder |
| stats | FOLDER | Created at dryrun but populated during run | Contains older timestamped `runinfo.yaml` files |
| submit_script.sbatch | TXT | run | Slurm script to kickstart the main Snakemake job |
| tools.yaml | YAML | run | YAML containing the version of tools used in the pipeline (obsolete; was used to load specific module versions prior to moving over to Docker/Singularity containers) |

## Resultsdir

The results directory contains the actual output files. Below are the folders that you may find within it.

```bash
WORKDIR
├── results
    ├── dedupBam
    ├── peaks
    ├── QC
    ├── qsortedBam
    ├── tagAlign
    └── tmp
```

Content details:

| Folder       | Description |
|-------------|------------|
| dedupBam    | Deduplicated filtered BAM files; can be used for visualization. |
| peaks       | Genrich/MACS2 peak calls (raw, consensus, fixed-width); also contains ROI files with Diff-ATAC results if `contrasts.tsv` is provided; motif enrichments using HOMER and AME; bigwigs for visualization. |
| QC          | Flagstats; dupmetrics; read counts; motif enrichments; FLD stats; Fqscreen; FRiP; ChIPSeeker results; TSS enrichments; Preseq; MultiQC. |
| qsortedBam  | Query name sorted BAM files; used for Genrich peak calling (includes multimappers). |
| tagAlign    | `tagAlign.gz` files; deduplicated; used for MACS2 peak calling. |
| tmp         | Can be deleted; blacklist index; intermediate FASTQs; Genrich output reads. |

Most of the above folders are self-explanatory. The `peaks` folder has this hierarchy:

```bash
WORKDIR
├── results
    ├── peaks
        ├── genrich
        │   ├── bigwig
        │   ├── <replicateName>.genrich.narrowPeak_motif_enrichment
        │   │   └── knownResults
        │   ├── DiffATAC
        │   ├── fixed_width
        │   └── tn5nicks
        └── macs2
            ├── bigwig
            ├── <replicateName>.macs2.narrowPeak_motif_enrichment
            │   └── knownResults
            ├── DiffATAC
            ├── fixed_width
            └── tn5nicks
```

`tn5nicks` folders host the per-replicate BAM files containing the Tn5 nicking sites in Genrich or MACS2 "peakcalling" reads, respectively. For easy visualization, they are converted to bigWig format and saved in respective `bigwig` folders. `DiffATAC` contains the DESeq2 differential accessiblity results, both per-contrast and aggregated accross all contrasts in `contrasts.tsv`. These results are solely based on tn5 nick counts.


> DISCLAIMER: This folder hierarchy is specific to v1.0.6 and is subject to change with version.