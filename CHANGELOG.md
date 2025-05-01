## ASPEN development version

## ASPEN 1.1.0

### General Changes

- Simplified the pipeline by removing unused or redundant rules, scripts, and parameters.
- Updated the pipeline to focus on core functionalities, including alignment, peak calling, differential accessibility, and quality control.

### Configuration Changes

- Added `spikein` and `spikein_genome` parameters to `config.yaml`.
- Updated `cluster.json` to increase the number of threads for most rules and removed unused rules.

### Alignment Changes

- Added spike-in alignment (`align2spikein`) and scaling factor computation rules.
- Updated `align.smk` to simplify alignment steps and add spike-in-related logic.
- Added the `chromosomes` parameter in alignment and replaced the hardcoded chromosome names.

### Peak Calling Changes

- Updated MACS2 and Genrich peak calling rules to add spike-in-related logic.
- Simplified peak filtering by applying a lenient q-value filter (default 0.1) directly in the peak calling scripts.
- Peak called narrowPeak files are now default q-value filtered. Unfiltered files are still available for users who want to apply their own filters.
- Removed redundant peak annotation steps and streamlined the output structure for peak calling.

### Differential Accessibility Changes

- Add new rules for scaling counts and annotating regions of interest.
- DiffATAC analysis is now run for both MACS2 and Genrich peak calls, with results stored in separate directories.
- Updated the `diffatac.smk` file to include spike-in scaling factors in the differential analysis.
- Removed redundant steps in the differential accessibility analysis to streamline the process.

### Quality Control Changes

- Updated FRiP calculation to use `tagAlign.gz` files instead of deduplicated BAM files.
- Removed unnecessary QC metrics and simplified the QC workflow.
- Updated TSS enrichment and fragment length distribution rules to align with the simplified pipeline structure.

### Documentation Changes

- Updated all documentation files to reflect the addition of spike-in normalization and related features.
- Simplified the documentation to focus on the core functionalities of the pipeline.
- Added references to spike-in genomes, scaling factors, and related analyses.

### Script Changes

- Appended scripts related to spike-in normalization, scaling factors, and Tn5 BAM generation.
- Updated peak calling scripts (`ccbr_atac_macs2_peak_calling.bash` and `ccbr_atac_genrich_peak_calling.bash`) to simplify logic and remove redundant steps.
- Updated featureCounts header fix script to remove references to spike-in-specific file naming.

### Output Structure Changes

- Simplified the output directory structure.
- Consolidated peak calling outputs into a single directory for each peak caller.
- Updated the `results` folder structure to align with the simplified pipeline.

### Miscellaneous Changes

- Updated all rules and scripts to use consistent naming conventions and file paths.
- Improved readability and maintainability of the codebase by removing redundant logic and comments.

## ASPEN 1.0.6

- fix: dockername typo ([#57](https://github.com/CCBR/ASPEN/issues/57), @kopardev)

## ASPEN 1.0.5

- fix: ucsc tool version changed requiring newer version of GLIBC ([#54](https://github.com/CCBR/ASPEN/issues/54), @kopardev)
- using new masterdocker v11

## ASPEN 1.0.4

- fix: DiffATAC failure ([#46](https://github.com/CCBR/ASPEN/issues/46), @kopardev)
- fix: last line of `contrasts.tsv` read in correctly; black lines ignored ([#48](https://github.com/CCBR/ASPEN/issues/48), @kopardev)
- fix: ROI calculation from fixed-width consensus peaks no longer tried to fix the peak width again ([#50](https://github.com/CCBR/ASPEN/issues/50), @kopardev)
- feature: create diffatac results from MACS2 peaks ([#51](https://github.com/CCBR/ASPEN/issues/51), @kopardev)
- fix: `BUYINPARTITIONS` fixed in wrapper for BIOWULF-only ([#52](https://github.com/CCBR/ASPEN/issues/52), @kopardev)

## ASPEN 1.0.3

- fix: No module named 'numpy.\_core.\_multiarray_umath' error with`unset PYTHONPATH` ([#43](https://github.com/CCBR/ASPEN/issues/43), @kopardev)
- fix: jobby command points to the correct location of snakemake.log file
- ASPEN is now archived on Zenodo, you can cite it with the DOI [10.5281/zenodo.13755867](https://doi.org/10.5281/zenodo.13755867). ([#42](https://github.com/CCBR/ASPEN/issues/42), @kelly-sovacool)

## ASPEN 1.0.2

- Set the singularity cache dir if `--singcache` is not provided. ([#37](https://github.com/CCBR/ASPEN/pull/37), @kelly-sovacool)
- ASPEN now has a documentation website: <https://ccbr.github.io/ASPEN>

## ASPEN 1.0.1

- differential ATAC updated
- documentation updated

## ASPEN 1.0.0

- completely dockerized
- differential ATAC

## ASPEN 0.6.1

- correction to fqscreen cattle path

## ASPEN 0.6

- support for mmul10 (Macaca) and bosTau9 (cattle) genomes
  - created resource files: indexes, promoter files, tss files etc.
- Added Macaca and Cattle to fastqscreen indexes
- support increased from 4 replicate to 6 replicates
- macs and genrich fixed width peaks generation rule added
- docker updated to v10 (genome support and tidyverse added)

## ASPEN 0.5.3

- Includes reference files for mmul10

## ASPEN 0.5.2

- atac_assign_multimappers.py now getting query sorted input
- dryrun log saved in workdir
- local (workdir) scriptsdir used

## ASPEN 0.5.1

- typo fix in main wrapper script

## ASPEN 0.5

- fastqscreen added
- minor bug fixes

## ASPEN 0.4.1

- Bug fixes
- minor updates

## ASPEN 0.4

- Multiqc edits
- README updates

## ASPEN 0.3

- FRiP calculations added

## ASPEN 0.2

- Peak motif enrichment with homer/meme
- Peak replicate/sample/peakcaller PCA comparisons after bedtools jaccard pairwise calculations

## ASPEN 0.1

- first working version
- calls peaks with macs2/genrich
- annotates peaks with chipseeker (human and mouse support)
