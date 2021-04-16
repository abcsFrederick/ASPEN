def get_fastqs(wildcards):
    d=dict()
    d["R1"]=SAMPLESDF["R1"][wildcards.sample]
    d["R2"]=SAMPLESDF["R2"][wildcards.sample]
    return d

rule trim_align_dedup:
    input:
        unpack(get_fastqs)
    output:
        ta=join(WORKDIR,"tagAlign","{sample}.tagAlign.gz"),
        fastqcraw1=join(WORKDIR,"QC","fastqc","{sample}.R1_fastqc.zip"),
        fastqcraw2=join(WORKDIR,"QC","fastqc","{sample}.R2_fastqc.zip"),
        fastqc1=join(WORKDIR,"QC","fastqc","{sample}.R1.noBL_fastqc.zip"),
        fastqc2=join(WORKDIR,"QC","fastqc","{sample}.R2.noBL_fastqc.zip"),
        nreads=join(WORKDIR,"QC","{sample}.nreads.txt"),
        nrf=join(WORKDIR,"QC","preseq","{sample}.nrf"),
        dedupbam=join(WORKDIR,"bam","{sample}.dedup.bam"),
        qsortedbam=join(WORKDIR,"bam","{sample}.qsorted.bam"),
        picardout=join(WORKDIR,"QC","{sample}.dupmetric"),
        genomefile=join(WORKDIR,"bam","{sample}.genome"),
        trimfq1=join(WORKDIR,"trim","{sample}.R1.trim.fastq.gz"),
        trimfq2=join(WORKDIR,"trim","{sample}.R2.trim.fastq.gz")
    params:
        genome=GENOME,
        indexdir=INDEXDIR,
        workdir=WORKDIR,
        scriptsdir=SCRIPTSDIR,
        qcdir=join(WORKDIR,"QC"),
        sample="{sample}",
        script="ccbr_atac_trim_align_pe.bash"
    container: config["masterdocker"]
    threads: 56
    shell:"""
ls {params.scriptsdir}/{params.script}
cd /lscratch/$SLURM_JOBID && bash {params.scriptsdir}/{params.script} \
--infastq1 {input.R1} \
--infastq2 {input.R2} \
--threads {threads} \
--genome {params.genome} \
--scriptsfolder {params.scriptsdir} \
--keepfiles True

ls -larth

rsync -az --progress {params.sample}.dedup.bam {params.workdir}/bam/
rsync -az --progress {params.sample}.genome {params.workdir}/bam/
rsync -az --progress {params.sample}.dedup.bam.bai {params.workdir}/bam/
rsync -az --progress {params.sample}.qsorted.bam {params.workdir}/bam/

rsync -az --progress {params.sample}.tagAlign.gz {params.workdir}/tagAlign/

rsync -az --progress {params.sample}.bowtie2.bam.flagstat {params.qcdir}/
rsync -az --progress {params.sample}.bowtie2.log {params.qcdir}/
rsync -az --progress {params.sample}.dedup.bam.flagstat {params.qcdir}/
rsync -az --progress {params.sample}.dupmetric {params.qcdir}/
rsync -az --progress {params.sample}.filt.bam.flagstat {params.qcdir}/
rsync -az --progress {params.sample}.nreads.txt {params.qcdir}/

rsync -az --progress {params.sample}.nrf {params.qcdir}/preseq/
rsync -az --progress {params.sample}.preseq {params.qcdir}/preseq/
rsync -az --progress {params.sample}.preseq.log {params.qcdir}/preseq/

rsync -az --progress {params.sample}.R1.noBL_fastqc.html {params.qcdir}/fastqc/
rsync -az --progress {params.sample}.R2.noBL_fastqc.html {params.qcdir}/fastqc/
rsync -az --progress {params.sample}.R1.noBL_fastqc.zip {params.qcdir}/fastqc/
rsync -az --progress {params.sample}.R2.noBL_fastqc.zip {params.qcdir}/fastqc/
rsync -az --progress {params.sample}.R1_fastqc.html {params.qcdir}/fastqc/
rsync -az --progress {params.sample}.R2_fastqc.html {params.qcdir}/fastqc/
rsync -az --progress {params.sample}.R1_fastqc.zip {params.qcdir}/fastqc/
rsync -az --progress {params.sample}.R2_fastqc.zip {params.qcdir}/fastqc/

rsync -az --progress {params.sample}.R1.trim.fastq.gz {params.workdir}/trim/
rsync -az --progress {params.sample}.R2.trim.fastq.gz {params.workdir}/trim/

"""


rule trim:
    # group: "TAD"
    input:
        unpack(get_fastqs)
    output:
        R1=join(WORKDIR,"tmp","trim","{sample}.R1.trim.fastq.gz"),
        R2=join(WORKDIR,"tmp","trim","{sample}.R2.trim.fastq.gz")
    params:
        sample="{sample}",
        scriptsdir=SCRIPTSDIR,
        script="ccbr_cutadapt_pe.bash"
    container: config["masterdocker"]    
    threads: getthreads("trim")
    shell:"""
cd /lscratch
bash {params.scriptsdir}/{params.script} \
--infastq1 {input.R1} \
--infastq2 {input.R2} \
--samplename {params.sample} \
--threads {threads} \
--outfastq1 {output.R1} \
--outfastq2 {output.R2}
"""

rule remove_BL:
    # group: "TAD"
    input:
        R1=rules.trim.output.R1,
        R2=rules.trim.output.R2
    output:
        R1=join(WORKDIR,"tmp","noBL","{sample}.R1.noBL.fastq.gz"),
        R2=join(WORKDIR,"tmp","noBL","{sample}.R2.noBL.fastq.gz")
    params:
        sample="{sample}",
        scriptsdir=SCRIPTSDIR,
        script="ccbr_remove_blacklisted_reads_pe.bash"
    container: config["masterdocker"]    
    threads: getthreads("remove_BL")
    shell:"""
cd /lscratch
bash {params.scriptsdir}/{params.script} \
--infastq1 {input.R1} \
--infastq2 {input.R2} \
--samplename {params.sample} \
--threads {threads} \
--outfastq1 {output.R1} \
--outfastq2 {output.R2}
"""