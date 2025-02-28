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
| `dedupBam`    | Deduplicated filtered BAM files; can be used for visualization. |
| `peaks`    | Genrich/MACS2 peak calls (raw, consensus, fixed-width); also contains ROI files with Diff-ATAC results if `contrasts.tsv` is provided; motif enrichments using HOMER and AME; bigwigs for visualization. |
| `QC`          | Flagstats; dupmetrics; read counts; motif enrichments; FLD stats; Fqscreen; FRiP; ChIPSeeker results; TSS enrichments; Preseq; MultiQC. |
| `qsortedBam`  | Query name sorted BAM files; used for Genrich peak calling (includes multimappers). |
| `tagAlign`    | `tagAlign.gz` files; deduplicated; used for MACS2 peak calling. |
| `tmp`         | Can be deleted; blacklist index; intermediate FASTQs; Genrich output reads. |

The `QC` folder contains the `multiqc_report.html` file which provides a comprehensive summary of the quality control metrics across all samples, including read quality, duplication rates, and other relevant statistics. This report aggregates results from various QC tools such as FastQC, FastqScreen, FLD, TSS enrichment, Peak Annotations, and others, presenting them in an easy-to-read format with interactive plots and tables. It helps in quickly identifying any issues with the sequencing data and ensures that the data quality is sufficient for downstream analysis.

!!! note
    BAM files from `dedupBam` can be used for downstream footprinting analysis using [CCBR_TOBIAS](https://github.com/CCBR/CCBR_Tobias) pipeline

!!! note
    [bamCompare](https://deeptools.readthedocs.io/en/develop/content/tools/bamCompare.html) from deeptools can be run to compare BAMs from `dedupBam` for comprehensive BAM comparisons.

!!! note
    BAM files from `dedupBam` can also be converted to BED format and processed with [chromVAR](https://github.com/GreenleafLab/chromVAR) to identify variability in motif accessibility across samples and assess differentially active transcription factors from the JASPAR database.

Most of the above folders are self-explanatory. The `peaks` folder has this hierarchy:

```bash
WORKDIR
├── results
    ├── peaks
        ├── genrich
        │   ├── <replicateName>.genrich.narrowPeak
        │   ├── <replicateName>.genrich.narrowPeak.annotated
        │   ├── <replicateName>.genrich.narrowPeak.genelist
        │   ├── <replicateName>.genrich.narrowPeak.annotation_summary
        │   ├── <replicateName>.genrich.narrowPeak.annotation_distribution
        │   ├── <sampleName>.genrich.pooled.narrowPeak
        │   ├── <sampleName>.genrich.consensus.bed
        │   ├── ROI.counts.tsv
        │   ├── bigwig
        │   ├── <replicateName>.genrich.narrowPeak_motif_enrichment
        │   │   └── knownResults
        │   ├── DiffATAC
        │   ├── <replicateName>.genrich.consensus.bed_motif_enrichment
        │   └── tn5nicks
        └── macs2
        │   ├── <replicateName>.macs2.narrowPeak
        │   ├── <replicateName>.macs2.narrowPeak.annotated
        │   ├── <replicateName>.macs2.narrowPeak.genelist
        │   ├── <replicateName>.macs2.narrowPeak.annotation_summary
        │   ├── <replicateName>.macs2.narrowPeak.annotation_distribution
        │   ├── <sampleName>.macs2.pooled.narrowPeak
        │   ├── <sampleName>.macs2.consensus.bed
        │   ├── ROI.counts.tsv
            ├── bigwig
            ├── <replicateName>.macs2.narrowPeak_motif_enrichment
            │   └── knownResults
            ├── DiffATAC
        │   ├── <replicateName>.macs2.consensus.bed_motif_enrichment
            ├── fixed_width
            └── tn5nicks
```

Some of the important folders and files are highlighted below:

### Folders

- `bigwig`: 

For easy visualization, they are converted to bigWig format and saved in respective `bigwig` folders. The bigWig files can be directly loaded into [UCSC Browser](https://genome.ucsc.edu/) or [IGV](https://igv.org/doc/desktop/).

- `tn5nicks`: 

This folder host the per-replicate BAM files containing the Tn5 nicking sites in Genrich or MACS2 "peakcalling" reads, respectively. 

- `DiffATAC`: 

Contains the DESeq2 differential accessiblity results, both per-contrast and aggregated accross all contrasts in `contrasts.tsv`. These results are solely based on tn5 nick counts.

- `fixed_width`: 

This folder contains fixed-width consensus peaks across replicates and samples, represented in the "Regions-Of-Interest" files. The `ROI.bed` file lists genomic regions where chromatin accessibility is analyzed using DESeq2, with results stored in the `DiffATAC` folder.

- `<replicateName>.macs2.narrowPeak_motif_enrichment`;`<replicateName>.genrich.narrowPeak_motif_enrichment`;`<replicateName>.macs2.consensus.bed_motif_enrichment`;`<replicateName>.genrich.consensus.bed_motif_enrichment`:

Contains the motif enrichments calculated using HOMER and AME for peaks called for each replicate, sample consensus peaks using both MACS2 and Genrich. Specifically, two types of motif enrichments are performed:
  
  - Enrichment of known [HOCOMOCO](https://hocomoco11.autosome.org/) (version 11) motifs for HUMAN or MOUSE or BOTH using [HOMER](http://homer.ucsd.edu/homer/ngs/peaks.html). See file `knownResults.html`.

  - _de novo_ motif enrichment using [AME](https://meme-suite.org/meme/doc/ame.html) from MEME suite. See file `ame_results.txt`. Custom parallelization is used to optimize AME based enrichment analysis.

### Files

- `*.narrowPeak`:

Called peaks from Genrich or MACS2

- Annotated peak files:

Peaks are annotated with ChIPSeeker and results are saved in the following files:

  - `.annotated`

    Tab-delimited txt file with the following columns:

| Column Number | Field Name | Description |
|--------------|------------|-------------|
| 1  | #peakID        | Peak identifier |
| 2  | chrom         | Peak chromosome |
| 3  | chromStart    | Peak start coordinate |
| 4  | chromEnd      | Peak end coordinate |
| 5  | width         | Peak width |
| 6  | annotation    | Peak annotation (Promoter; 3' or 5' UTR; Distal; Downstream; Exon; Intron) |
| 7  | geneChr       | Gene chromosome |
| 8  | geneStart     | Gene start coordinate |
| 9  | geneEnd       | Gene end coordinate |
| 10 | geneLength    | Gene length (including introns) |
| 11 | geneStrand    | Gene strand |
| 12 | geneId        | Gene identifier |
| 13 | transcriptId  | Transcript identifier |
| 14 | distanceToTSS | Distance of peak from the Transcription Start Site |
| 15 | ENSEMBL       | Gene Ensembl ID |
| 16 | SYMBOL        | Gene symbol |
| 17 | GENENAME      | Gene description |
| 18 | score         | Score from `.narrowPeak` file |
| 19 | signalValue   | Signal from `.narrowPeak` file |
| 20 | pValue        | p-value from `.narrowPeak` file |
| 21 | qValue        | q-value from `.narrowPeak` file |
| 22 | peak          | Distance of peak summit from peak start coordinate |

  - `.genelist`

This is a tab-delimited file with names (Ensembl ID, gene symbol) of genes which have ATAC-seq peaks in their promotor regions. This file can be used downstream for gene enrichment analysis (ORA or over-representation analysis).

  - `.annotation_summary`; `.annotation_distribution`

Tab-delimited files that provide statistics on peak annotations, quantifying the number of peaks found in Promoters, Exonic regions, Distal Intergenic regions, etc. The `.annotation_distribution` is use to create visualization of these annotation-distributions in the MultiQC report.

- `ROI.counts.tsv`

This file contains the read counts for each Region-Of-Interest (ROI) across all replicates of all samples. It is a tab-delimited file with the following columns:

| Column Number | Field Name | Description |
|---------------|------------|-------------|
| 1             | Geneid       | Region-Of-Interest identifier |
| 2             | Chr      | Chromosome of the ROI |
| 3             | Start      | Start coordinate of the ROI |
| 4             | End        | End coordinate of the ROI |
| 5             | Strand        | "." |
| 6             | Length     | Length of the ROI |
| 7             | sample1_replicate1    | Tn5 nicking site counts in this ROI for replicate1 of sample1 |
| 8             | sample1_replicate2    | Tn5 nicking site counts in this ROI for replicate2 of sample1 |
| ...           | ...        | ... |
| n             | sampleN_replicateM    | Tn5 nicking site counts in this ROI for replicateM of sampleN|

Each row represents a specific ROI, and the columns contain the read counts for each sample, allowing for differential accessibility analysis.

!!! warning
    DISCLAIMER: This folder hierarchy is specific to v1.0.6 and is subject to change with version.