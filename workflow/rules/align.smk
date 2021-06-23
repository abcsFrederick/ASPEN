def get_fastqs(wildcards):
    d=dict()
    d["R1"]=REPLICATESDF["R1"][wildcards.replicate]
    d["R2"]=REPLICATESDF["R2"][wildcards.replicate]
    return d

localrules: align_stats


rule trim:
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

rule remove_BL:
    # group: "TAD"
    input:
        R1=rules.trim.output.R1,
        R2=rules.trim.output.R2,
    output:
        R1=join(RESULTSDIR,"tmp","trim","{replicate}.R1.noBL.fastq.gz"),
        R2=join(RESULTSDIR,"tmp","trim","{replicate}.R2.noBL.fastq.gz"),
    params:
        replicate="{replicate}",
        scriptsdir=SCRIPTSDIR,
        genome=GENOME,
        script="ccbr_remove_blacklisted_reads_pe.bash"
    container: config["masterdocker"]    
    threads: getthreads("remove_BL")
    shell:"""
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi
bash {params.scriptsdir}/{params.script} \
--infastq1 {input.R1} \
--infastq2 {input.R2} \
--samplename {params.replicate} \
--genome {params.genome} \
--threads {threads} \
--outfastq1 {output.R1} \
--outfastq2 {output.R2}
"""

rule align:
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
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi
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

rule align_stats:
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
