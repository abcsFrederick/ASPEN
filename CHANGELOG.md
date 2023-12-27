## v1.0.1

-   differential ATAC updated
-   documentation updated

## v1.0.0

-   completely dockerized
-   differential ATAC

## v0.6.1 Latest

-   correction to fqscreen cattle path

## v0.6

-   support for mmul10 (Macaca) and bosTau9 (cattle) genomes
    -   created resource files: indexes, promoter files, tss files etc.
-   Added Macaca and Cattle to fastqscreen indexes
-   support increased from 4 replicate to 6 replicates
-   macs and genrich fixed width peaks generation rule added
-   docker updated to v10 (genome support and tidyverse added)

## v0.5.3

-   Includes reference files for mmul10

## v0.5.2

-   atac_assign_multimappers.py now getting query sorted input
-   dryrun log saved in workdir
-   local (workdir) scriptsdir used

## v0.5.1

-   typo fix in main wrapper script

## v0.5

-   fastqscreen added
-   minor bug fixes

## v0.4.1

-   Bug fixes
-   minor updates

## v0.4

-   Multiqc edits
-   README updates

## v0.3

-   FRiP calculations added

## v0.2

-   Peak motif enrichment with homer/meme
-   Peak replicate/sample/peakcaller PCA comparisons after bedtools jaccard pairwise calculations

## v0.1

-   first working version
-   calls peaks with macs2/genrich
-   annotates peaks with chipseeker (human and mouse support)
