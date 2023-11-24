#########################################################
# FUNCTION DEFINITIONS
#########################################################

def get_fastqs(wildcards):
    """
    Return a dictionary of input files for the trimming rule.
    """
    d=dict()
    d["R1"]=REPLICATESDF["R1"][wildcards.replicate]
    d["R2"]=REPLICATESDF["R2"][wildcards.replicate]
    return d

#########################################################

localrules: align_stats

#########################################################

rule trim:
# """
# First step of TAD
# Input: raw fastq files
# Output: trimmed fastq files in "tmp" folder, so they can be
# deleted if no longer required.
# """
    # group: "TAD"
    input:
        unpack(get_fastqs)
    output:
        R1=join(RESULTSDIR,"tmp","trim","{replicate}.R1.trim.fastq.gz"),
        R2=join(RESULTSDIR,"tmp","trim","{replicate}.R2.trim.fastq.gz"),
    params:
        replicate="{replicate}",
        scriptsdir=SCRIPTSDIR,
        script="ccbr_cutadapt_pe.bash"
    container: config["masterdocker"]
    threads: getthreads("trim")
    shell:"""
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi
bash {params.scriptsdir}/{params.script} \
--infastq1 {input.R1} \
--infastq2 {input.R2} \
--samplename {params.replicate} \
--threads {threads} \
--outfastq1 {output.R1} \
--outfastq2 {output.R2}
"""

#########################################################

rule create_BL_index:
# """
# Create bwa index from blacklist fasta file
# """
    # group: "TAD"
    input:
        fa=BLACKLISTFA,
    output:
        sa=join(RESULTSDIR,"tmp","BL",GENOME+"_blacklist.sa")
    params:
        # bldir=join(RESULTSDIR,"tmp","BL"),
        genome=GENOME,
    container: config["masterdocker"]
    threads: getthreads("create_BL_index")
    shell:"""
blindexdir=$(dirname {output.sa})
sabasename=$(basename {output.sa})
indexprefix=${{sabasename%.*}}
cd $blindexdir
bwa index -p $indexprefix {input.fa}
"""

#########################################################

rule remove_BL:
# """
# Second step of TAD
# Input: trimmed fastq files
# Output: not aligning to blacklisted regions, trimmed fastq files
# in "tmp" folder, so they can be deleted if no longer required.
# """
    # group: "TAD"
    input:
        R1=rules.trim.output.R1,
        R2=rules.trim.output.R2,
        sa=rules.create_BL_index.output.sa
    output:
        R1=join(RESULTSDIR,"tmp","trim","{replicate}.R1.noBL.fastq.gz"),
        R2=join(RESULTSDIR,"tmp","trim","{replicate}.R2.noBL.fastq.gz"),
    params:
        replicate="{replicate}",
        scriptsdir=SCRIPTSDIR,
        script="ccbr_remove_blacklisted_reads_pe.bash"
    container: config["masterdocker"]
    threads: getthreads("remove_BL")
    shell:"""
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi
sa={input.sa}
blindex=${{sa%.*}}
bash {params.scriptsdir}/{params.script} \
--infastq1 {input.R1} \
--infastq2 {input.R2} \
--samplename {params.replicate} \
--blindexdir $blindex \
--threads {threads} \
--outfastq1 {output.R1} \
--outfastq2 {output.R2}
"""

#########################################################

rule align:
# """
# Third (final) step of TAD
# Alignment is handled by bowtie
# Filtering is based off of ENCODEs pipeline (atac_assign_multimappers.py is directly pulled from ENCODE)
# Input: not aligning to blacklisted regions, trimmed fastq files
# Output:
#     1. aligned and filtered alignment files:
#         * tagAlign --> for MACS2
#         * dedup bam --> for TOBIAS
#         * qsorted bam --> for Genrich
#     2. QC files:
#         * nrf --> non-redundant fraction using preseq
#         * dupmetric --> to quantify duplication level using picard
#         * flagstat files --> quatify surviving reads at each step
# """
    # group: "TAD"
    input:
        R1=rules.remove_BL.output.R1,
        R2=rules.remove_BL.output.R2,
    output:
        tagAlign=join(RESULTSDIR,"tagAlign","{replicate}.tagAlign.gz"),
        dedupbam=join(RESULTSDIR,"dedupBam","{replicate}.dedup.bam"),
        qsortedBam=join(RESULTSDIR,"qsortedBam","{replicate}.qsorted.bam"),
        fs1=join(QCDIR,"{replicate}.bowtie2.bam.flagstat"),
        fs2=join(QCDIR,"{replicate}.dedup.bam.flagstat"),
        fs3=join(QCDIR,"{replicate}.filt.bam.flagstat"),
        dupmetric=join(QCDIR,"{replicate}.dupmetric"),
        nrf=join(QCDIR,"preseq","{replicate}.nrf"),
    params:
        replicate="{replicate}",
        workdir=RESULTSDIR,
        qcdir=QCDIR,
        indexdir=INDEXDIR,
        scriptsdir=SCRIPTSDIR,
        genome=GENOME,
        script="ccbr_bowtie2_align_pe.bash",
        multimapping=config["multimapping"],
        mem=getmemG("align")
    container: config["masterdocker"]
    threads: getthreads("align")
    shell:"""
set -exo pipefail

if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi

# make folder if they do not exist
for  f in {output.tagAlign} {output.dedupbam} {output.qsortedBam};do
if [ ! -w $(dirname $f) ];then
mkdir -p $(dirname $f)
fi
done

bash {params.scriptsdir}/{params.script} \
--infastq1 {input.R1} \
--infastq2 {input.R2} \
--samplename {params.replicate} \
--genomename {params.genome} \
--threads {threads} \
--mem {params.mem} \
--indexdir {params.indexdir} \
--multimapping {params.multimapping} \
--scriptsfolder {params.scriptsdir}

rsync -az --progress --verbose --remove-source-files {params.replicate}.tagAlign.gz {output.tagAlign}
rsync -az --progress --verbose --remove-source-files {params.replicate}.qsorted.bam {output.qsortedBam}


rsync -az --progress --verbose --remove-source-files {params.replicate}.bowtie2.bam.flagstat {params.qcdir}/
rsync -az --progress --verbose --remove-source-files {params.replicate}.bowtie2.log {params.qcdir}/
rsync -az --progress --verbose --remove-source-files {params.replicate}.dedup.bam {output.dedupbam}
rsync -az --progress --verbose --remove-source-files {params.replicate}.dedup.bam.bai $(dirname {output.dedupbam})/

rsync -az --progress --verbose --remove-source-files {params.replicate}.dedup.bam.flagstat {params.qcdir}/
rsync -az --progress --verbose --remove-source-files {params.replicate}.dupmetric {params.qcdir}/
rsync -az --progress --verbose --remove-source-files {params.replicate}.filt.bam.flagstat {params.qcdir}/

rsync -az --progress --verbose --remove-source-files {params.replicate}.nrf {params.qcdir}/preseq/
rsync -az --progress --verbose --remove-source-files {params.replicate}.preseq {params.qcdir}/preseq/
rsync -az --progress --verbose --remove-source-files {params.replicate}.preseq.log {params.qcdir}/preseq/
"""

#########################################################

rule align_stats:
# """
# Collect stats at each TAD step
# Output: per-replicate nreads.txt file --> required by MultiQC for reporting.
# """
    # group: "TAD"
    input:
        unpack(get_fastqs),
        trimR1=rules.trim.output.R1,
        noBLR1=rules.remove_BL.output.R1,
        fs1=rules.align.output.fs3,
        fs2=rules.align.output.fs2
    output:
        nreads=join(QCDIR,"{replicate}.nreads.txt"),
    params:
        replicate="{replicate}",
    shell:"""
nreads=`zcat {input.R1}|wc -l`
nreadstrim=`zcat {input.trimR1}|wc -l`
echo "$nreads $nreadstrim"|awk '{{printf("%d\\tInput Nreads\\n%d\\tAfter trimming\\n",$1/2,$2/2)}}' > {output.nreads}
nreads=`zcat {input.noBLR1}|wc -l`
echo "$nreads"|awk '{{printf("%d\\tAfter removing mito-ribo reads\\n",$1/2)}}' >> {output.nreads}
nreads=`grep -m1 total {input.fs1}|awk '{{print $1}}'`
echo "$nreads"|awk '{{printf("%d\\tMapped reads\\n",$1)}}' >> {output.nreads}
nreads=`grep -m1 total {input.fs2}|awk '{{print $1}}'`
echo "$nreads"|awk '{{printf("%d\\tAfter deduplication\\n",$1)}}' >> {output.nreads}
"""

#########################################################
