QCDIR=join(RESULTSDIR,"QC")


rule fastqc:
    input:
        expand(join(WORKDIR,"fastqs","{sample}.R1.fastq.gz"),sample=SAMPLES),
        expand(join(WORKDIR,"fastqs","{sample}.R2.fastq.gz"),sample=SAMPLES),
        expand(rules.remove_BL.output.R1,sample=SAMPLES),
        expand(rules.remove_BL.output.R2,sample=SAMPLES),
    output:
        expand(join(QCDIR,"fastqc","{sample}.R1_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R2_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R1.noBL_fastqc.zip"), sample=SAMPLES),
        expand(join(QCDIR,"fastqc","{sample}.R2.noBL_fastqc.zip"), sample=SAMPLES),
    params:
        outdir=join(QCDIR,"fastqc"),
    threads: getthreads("fastqc")
    container: config["masterdocker"]
    shell: """
    fastqc {input} -t {threads} -o {params.outdir};
    """      


rule atac_tss:
    input:
        tagalign=join(RESULTSDIR,"tagAlign","{sample}.tagAlign.gz")
    output:
        tss=join(QCDIR,"tss","{sample}.tss.txt")
    params:
        sample="{sample}",
        workdir=RESULTSDIR,
        qcdir=join(QCDIR),
        indexdir=INDEXDIR,
        scriptsdir=SCRIPTSDIR,
        tssdir=join(RESOURCESDIR,"db"),
        genome=GENOME,
        script="ccbr_tagAlign2TSS.bash",
    container: config["masterdocker"]    
    threads: getthreads("align")
    shell:"""
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi
tagAlign=$(basename {input.tagalign})
outfn=$(basename {output.tss})

rsync -az --progress --remove-source-files {input.tagalign} $tagAlign

bash {params.scriptsdir}/{params.script} \
--tagaligngz $tagAlign \
--genome {params.genome} \
--tsstxt $outfn \
--ncpus {threads} \
--scriptsfolder {params.scriptsdir} \
--resourcesfolder {params.tssdir}

rsync -az --progress --remove-source-files $outfn {output.tss} && rm -f *
"""
