#!/bin/bash

# . /opt2/conda/etc/profile.d/conda.sh
# conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Adapter trimmed (and blacklist filtered) fastqs are aligned to genome using bowtie2, multimappers are properly assigned, deduplicated using picard, filtered based on mapq, bams converted to tagAlign files."      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--samplename',required=True, help='samplename')
parser.add_argument('--threads',required=True, help='number of threads')
parser.add_argument('--mem',required=True, help='mem in G for java')
parser.add_argument('--genomename',required=True, help='hg19/hg38/mm9/mm10')
parser.add_argument('--indexdir',required=True, help='dir where the bowtie2 index resides')
parser.add_argument('--multimapping',required=False, default=4, help='multimapping threshold')
parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
parser.add_argument('--keepfiles',required=False,default='False',help='to keep intermediate files, set this to True')
EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
samplename=$SAMPLENAME
genome=$GENOMENAME
multimapping=$MULTIMAPPING
log="${samplename}.bowtie2.log"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
ncpus=$THREADS

bowtie2 -X2000 -k $multimapping --very-sensitive --threads $ncpus -x ${INDEXDIR}/$genome \
 -1 $infastq1 -2 $infastq2 > ${samplename}.bowtie2.sam \
 2> $log


samtools view -@ $ncpus -bS -o ${samplename}.bowtie2.bam ${samplename}.bowtie2.sam
rm -rf ${samplename}.bowtie2.sam

samtools sort -@ $ncpus -o ${samplename}.bowtie2.sorted.bam ${samplename}.bowtie2.bam
mv ${samplename}.bowtie2.sorted.bam ${samplename}.bowtie2.bam
samtools flagstat ${samplename}.bowtie2.bam > ${samplename}.bowtie2.bam.flagstat
samtools index ${samplename}.bowtie2.bam

# -F 516
# remove reads that
# 1. are unmapped
# 2. fail platform/vendor qualit checks
samtools view -@ $ncpus -F 516 -u ${samplename}.bowtie2.bam $CHROMOSOMES > ${samplename}.tmp1.bam
samtools sort -@ $ncpus -n -o ${samplename}.tmp1.sorted.bam ${samplename}.tmp1.bam
mv ${samplename}.tmp1.sorted.bam ${samplename}.qsorted.bam
rm -rf ${samplename}.tmp1.bam
# qsorted.bam is used for peak calling with Genrich

# assign multimappers using ENCODE script
samtools view -@ $ncpus -h ${samplename}.qsorted.bam | \
${SCRIPTSFOLDER}/atac_assign_multimappers.py -k $multimapping --paired-end | \
samtools view -@ $ncpus -bS -o ${samplename}.tmp3.bam -

samtools sort -@ $ncpus -o ${samplename}.tmp3.sorted.bam ${samplename}.tmp3.bam
mv ${samplename}.tmp3.sorted.bam ${samplename}.tmp3.bam

# calculate non-redundant fraction before deduplication using PRESEQ
bash ${SCRIPTSFOLDER}/ccbr_bam2nrf.bash \
--bam ${samplename}.tmp3.bam \
--preseq ${samplename}.preseq \
--preseqlog ${samplename}.preseq.log \
--nrf ${samplename}.nrf \
--scriptsfolder $SCRIPTSFOLDER

# -F 256
# remove all non-primary alignments
samtools view -@ $ncpus -F 256 -u ${samplename}.tmp3.bam > ${samplename}.tmp4.bam
samtools sort -@ $ncpus -o ${samplename}.dup.bam ${samplename}.tmp4.bam
rm -rf ${samplename}.tmp3.bam ${samplename}.tmp4.bam

# -F 1796
# remove reads that are:
# 1. unmapped
# 2. not primary alignments
# 3. fail platform/vendor quality checks
# 4. is a PCR or optical duplicate
samtools view -@ $ncpus -F 1796 -u ${samplename}.dup.bam > ${samplename}.tmp5.bam
samtools sort -@ $ncpus -o ${samplename}.filt.bam ${samplename}.tmp5.bam
rm -rf ${samplename}.dup.bam ${samplename}.tmp5.bam

samtools index ${samplename}.filt.bam
samtools flagstat ${samplename}.filt.bam > ${samplename}.filt.bam.flagstat

# estimate duplication rate using picard
java -Xmx${MEM} -jar /opt2/picardcloud.jar MarkDuplicates \
INPUT=${samplename}.filt.bam \
OUTPUT=${samplename}.dupmark.bam \
METRICS_FILE=${samplename}.dupmetric \
VALIDATION_STRINGENCY=LENIENT \
ASSUME_SORTED=true \
REMOVE_DUPLICATES=false

ls -alrth

if [ $KEEPFILES == "False" ];then rm -rf ${samplename}.filt.bam;fi

# remove duplicates marked by picard to create dedup.bam from dupmark.bam
samtools view -F 1796 -b -@ $ncpus -o ${samplename}.dedup.tmp.bam ${samplename}.dupmark.bam
if [ $KEEPFILES == "False" ];then rm -rf ${samplename}.dupmark.bam;fi

samtools index ${samplename}.dedup.tmp.bam

# . /opt2/conda/etc/profile.d/conda.sh
# filter all alignments under MAPQ of 6 (5 or lower)
python ${SCRIPTSFOLDER}/ccbr_bam_filter_by_mapq.py -i ${samplename}.dedup.tmp.bam -o ${samplename}.dedup.bam -q 6
rm -rf ${samplename}.dedup.tmp.bam

samtools index ${samplename}.dedup.bam
samtools flagstat ${samplename}.dedup.bam > ${samplename}.dedup.bam.flagstat
samtools view -H ${samplename}.dedup.bam|grep "^@SQ"|cut -f2,3|sed "s/SN://g"|sed "s/LN://g" > ${samplename}.genome

# apply Tn5 shifts and create tagAlign files for peak calling with MACS2
bedtools bamtobed -i ${samplename}.dedup.bam | \
awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'| \
awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | \
gzip -c > ${samplename}.tagAlign.gz



# conda deactivate
