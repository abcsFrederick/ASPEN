# Overview of ASPEN

To address [these](introduction.md#challenges-in-atac-seq-data-analysis) challenges, CCBR developed ASPEN (**A**tac **S**eq **P**ip**E**li**N**e), an automated pipeline tailored for the comprehensive analysis of ATAC-seq data. ASPEN is designed to process paired-end Illumina sequencing data, guiding users from raw data through quality control, alignment, peak calling, and downstream analyses. By integrating a suite of established bioinformatics tools and adhering to best practices, ASPEN ensures robust, reproducible, and high-quality results. For more information on the challenges in ATAC-seq data analysis, refer to the [Challenges in ATAC-seq Data Analysis] section.

## **Data Preprocessing**

### **Adapter Trimming** 

Utilizes CutAdapt to remove adapter sequences from raw reads, ensuring that subsequent analyses are not confounded by extraneous sequences.

### **Quality Assessment**

Employs FastQC for initial quality checks and MultiQC for aggregated reporting, providing comprehensive insights into data quality and highlighting potential issues such as low-quality bases or GC content biases.


## **Alignment and Filtering**

### **Blacklist Filtering** 

Removes reads aligning to known blacklisted regions and mitochondria, reducing potential biases in peak detection and ensuring that analyses focus on biologically relevant regions.

### **Read Alignment** 

Aligns reads to the reference genome using Bowtie2, generating multiple output files, including `qsorted.bam`/`tagAlign.gz` for peak calling and `dedup.bam` for downstream analyses.


## **Peak Calling and Annotation**

### **Peak Detection** 

Implements MACS2 and Genrich to identify regions of open chromatin, accommodating both narrow and broad peak detection to capture a wide range of regulatory elements. Default MACS2 are selected from ENCODE. Peak calling uses lenient q-value setting and a q-value filter is then applied on the `.narrowPeak` file.

### **Consensus Peaks** 

For datasets with multiple replicates, ASPEN generates consensus peak sets based on simple overlaps between replicates and keeping peaks represented in >= 50% of replicates. Thus enhancing the reliability and reproducibility of detected regions by focusing on peaks consistently observed across replicates.

### **Peak Annotation** 

Utilizes ChIPseeker to annotate peaks, providing context regarding their genomic locations and potential regulatory roles, such as associations with promoters, enhancers, or gene bodies.

## **Quality Control Metrics**

### **Fragment Length Distribution**

Custom scripts are utilized to analyze insert size to assess nucleosome positioning and sample quality, with characteristic patterns indicating the presence of nucleosome-free regions and mono- or di-nucleosomes.

### **Library Complexity**

Estimates complexity using Preseq, aiding in the evaluation of sequencing depth sufficiency and potential PCR biases.

### **Transcription Start Site (TSS) Enrichment**

Calculates TSS enrichment scores, serving as indicators of data quality by measuring the accumulation of reads around transcription start sites, a hallmark of open chromatin.

### **Fraction of Reads in Peaks (FRiP)**

Determines the proportion of reads within called peaks, reflecting the signal-to-noise ratio and overall quality of the dataset. Fractions if reads in Promoters, Enhancers and DHS (DNase Hyper Sensitivity) regions are also computed.

## **Motif Enrichment Analysis**

### **HOMER and AME** 

Conducts motif enrichment analyses to identify potential transcription factor binding sites within accessible regions, offering insights into regulatory mechanisms and facilitating the discovery of key drivers of gene expression.

## **Comprehensive Reporting**

### **MultiQC Integration** 

Generates an extensive HTML report consolidating all QC metrics and analysis results, facilitating easy interpretation, visualization, and sharing of findings.

## **Differential Accessiblity Analyses**

In addition to its core functionalities, ASPEN offers advanced capabilities for differential accessibility analysis, enabling researchers to identify and interpret regions of the genome that exhibit significant changes in chromatin accessibility under different conditions. This process involves several key steps: refining peak calls to fixed-width peaks, defining regions of interest (ROIs), quantifying Tn5 transposase insertion event counts, performing statistical analysis to detect differential signals, and integrating these findings with ChIP-seq annotations for comprehensive biological insights.

### **Refining Peak Calls to Fixed-Width Peaks**

After the initial peak calling, ASPEN refines these peaks by centering them on their summits and extending them to a uniform width of 500 base pairs (bp). This standardization facilitates consistent downstream analyses and comparisons across different datasets. The approach of extending peak summits by ±250 bp to achieve a fixed width of 500 bp is a common practice in ATAC-seq data analysis, as it focuses on the most accessible regions of the chromatin and reduces variability in peak sizes. 


### **Defining Regions of Interest (ROIs)**

ASPEN processes fixed width peak calls on each replicate across all samples. Subsequently, it generates a union of all identified peaks to create comprehensive regions of interest (ROIs). The pipeline produces `ROI.bed` and `ROI.gtf` files for each of the peak callers, which are utilized for quantifying counts using featureCounts. This approach ensures a consistent framework for downstream differential accessibility analysis, enhancing the robustness and interpretability of the results.

### **Quantification of Differential Chromatin Accessibility**

Once the ROIs are established, ASPEN generates Tn5 nicking sites count matrices:

#### **Tn5 Nicking Sites Count Matrix** 

This matrix quantifies the frequency of Tn5 transposase insertion events at each ROI. Reads used by MACS2 and reads coming out of Genrich are used to compute Tn5 sites and counting them. The Tn5 transposase preferentially inserts into accessible regions of the chromatin, and the number of insertion events serves as a proxy for chromatin accessibility. By counting these insertion sites, researchers can accurately infer the openness of chromatin regions under different experimental conditions.


#### **Differential Accessibility Analysis with DESeq2**

To identify statistically significant differences in chromatin accessibility between conditions, ASPEN employs DESeq2, a widely-used tool for differential analysis of count data. DESeq2 models the count data to detect changes in accessibility, providing robust statistical inferences even with complex experimental designs. By comparing the count matrices across conditions, DESeq2 identifies ROIs with significant differential accessibility, shedding light on regulatory elements that may be functionally relevant to the experimental context.

#### **Integration with Gene Annotations**

To enhance the biological interpretation of differential accessibility results, ASPEN integrates the findings with peak and gene annotations from ChIPSeeker, offering insights into the regulatory landscape of the genome. This integration allows researchers to identify genes located near differentially accessible peaks, which can be crucial for understanding the functional implications of chromatin accessibility changes. By associating these peaks with nearby genes, scientists can infer potential regulatory relationships and gain a deeper understanding of the underlying molecular mechanisms.

## **Reporting**

ASPEN enhances differential chromatin accessibility analysis by providing an interactive HTML report generated from DESeq2 results, featuring various visualizations. Additionally, it offers a TSV (tabl-delimited) file with integrated gene annotations, compatible with Microsoft Excel, enabling efficient data manipulation and facilitating the identification of genes near regions with altered accessibility. The Excel file is aggregated accross all different contrasts queried in the project.

