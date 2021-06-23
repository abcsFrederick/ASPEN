def get_replicate_tagAlignFiles_for_sample(wildcards):
    replicates=SAMPLE2REPLICATES[wildcards.sample]
    d=dict()
    for i,s in enumerate(replicates):
        d["tagAlign"+str(i+1)]=join(RESULTSDIR,"tagAlign",s+".tagAlign.gz")
    return d

def get_replicate_qsortedBamFiles_for_sample(wildcards):
    replicates=SAMPLE2REPLICATES[wildcards.sample]
    d=dict()
    for i,s in enumerate(replicates):
        d["qsortedBam"+str(i+1)]=join(RESULTSDIR,"qsortedBam",s+".qsorted.bam")
    return d

rule atac_macs_peakcalling:
    input:
        unpack(get_replicate_tagAlignFiles_for_sample)
    params:
        genome=GENOME,
        genomefile=GENOMEFILE,
        outdir=join(RESULTSDIR,"peaks","macs2"),
        scriptsdir=SCRIPTSDIR,
        qcdir=QCDIR,
        script="ccbr_atac_macs2_peak_calling.bash",
        sample="{sample}",
        macs2_extsize=config["macs2"]["extsize"],
        macs2_shiftsize=config["macs2"]["shiftsize"]
    output:
        consensusPeakFileList=join(RESULTSDIR,"peaks","macs2","{sample}.consensus.macs2.peakfiles"),
        replicatePeakFileList=join(RESULTSDIR,"peaks","macs2","{sample}.replicate.macs2.peakfiles"),
        tn5nicksFileList=join(RESULTSDIR,"peaks","macs2","{sample}.macs2.tn5nicksbedfiles")
    container: config["masterdocker"]    
    threads: getthreads("trim")
    shell:"""
set -e -x -o pipefail
nreplicates=0
for f in {input}
do 
    nreplicates=$((nreplicates+1))
done

cd {params.outdir}

if [ ! -d {params.outdir}/bigwig ];then mkdir -p {params.outdir}/bigwig ;fi
if [ ! -d {params.outdir}/tn5nicks ];then mkdir -p {params.outdir}/tn5nicks ;fi
if [ ! -d {params.qcdir}/peak_annotation ];then mkdir -p {params.qcdir}/peak_annotation;fi

bash {params.scriptsdir}/{params.script} \
    --tagalignfiles {input} \
    --genome {params.genome} \
    --genomefile {params.genomefile} \
    --samplename {params.sample} \
    --scriptsfolder {params.scriptsdir} \
    --outdir {params.outdir} \
    --extsize {params.macs2_extsize} \
    --shiftsize {params.macs2_shiftsize}


for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
do
    for f in `ls *.${{g}}`;do
        rsync -az --progress  --remove-source-files $f {params.qcdir}/peak_annotation/
    done
done

if [ -f {output.tn5nicksFileList} ];then rm -f {output.tn5nicksFileList};fi
if [ -f {output.replicatePeakFileList} ];then rm -f {output.replicatePeakFileList};fi

for f in {input}
do
    repname=`basename $f|awk -F".tagAlign" '{{print $1}}'`
    echo -ne "${{repname}}\t{params.sample}\t{params.outdir}/tn5nicks/${{repname}}.macs2.tn5nicks.bam\n" >> {output.tn5nicksFileList}
    echo -ne "${{repname}}\t{params.sample}\t{params.outdir}/${{repname}}.macs2.narrowPeak\n" >> {output.replicatePeakFileList}
done
echo -ne "{params.sample}\t{params.sample}\t{params.outdir}/{params.sample}.macs2.consensus.bed\n" > {output.consensusPeakFileList}

"""


rule atac_genrich_peakcalling:
    input:
        unpack(get_replicate_qsortedBamFiles_for_sample)
    params:
        genome=GENOME,
        genomefile=GENOMEFILE,
        outdir=join(RESULTSDIR,"peaks","genrich"),
        scriptsdir=SCRIPTSDIR,
        qcdir=QCDIR,
        script="ccbr_atac_genrich_peak_calling.bash",
        sample="{sample}",
        genrich_s=config["genrich"]["s"],
        genrich_m=config["genrich"]["m"],
        genrich_q=config["genrich"]["q"],
        genrich_l=config["genrich"]["l"],
        genrich_g=config["genrich"]["g"],
    output:
        consensusPeakFileList=join(RESULTSDIR,"peaks","genrich","{sample}.consensus.genrich.peakfiles"),
        replicatePeakFileList=join(RESULTSDIR,"peaks","genrich","{sample}.replicate.genrich.peakfiles"),
        tn5nicksFileList=join(RESULTSDIR,"peaks","genrich","{sample}.genrich.tn5nicksbedfiles")
    container: config["masterdocker"]    
    threads: getthreads("trim")
    shell:"""
set -e -x -o pipefail
nreplicates=0
for f in {input}
do 
    nreplicates=$((nreplicates+1))
done

cd {params.outdir}

if [ ! -d {params.outdir}/bigwig ];then mkdir -p {params.outdir}/bigwig ;fi
if [ ! -d {params.outdir}/tn5nicks ];then mkdir -p {params.outdir}/tn5nicks ;fi
if [ ! -d {params.qcdir}/peak_annotation ];then mkdir -p {params.qcdir}/peak_annotation;fi

bash {params.scriptsdir}/{params.script} \
    --bamfiles {input} \
    --genome {params.genome} \
    --genomefile {params.genomefile} \
    --samplename {params.sample} \
    --scriptsfolder {params.scriptsdir} \
    --genrich_s {params.genrich_s} \
    --genrich_m {params.genrich_m} \
    --genrich_q {params.genrich_q} \
    --genrich_l {params.genrich_l} \
    --genrich_g {params.genrich_g} \
    --outdir {params.outdir}

for g in "annotated" "genelist" "annotation_summary" "annotation_distribution"
do
    for f in `ls *.${{g}}`;do
        rsync -az --progress  --remove-source-files $f {params.qcdir}/peak_annotation/
    done
done

if [ -f {output.tn5nicksFileList} ];then rm -f {output.tn5nicksFileList};fi
if [ -f {output.replicatePeakFileList} ];then rm -f {output.replicatePeakFileList};fi

for f in {input}
do
    repname=`basename $f|awk -F".bam" '{{print $1}}'`
    echo -ne "${{repname}}\t{params.sample}\t{params.outdir}/tn5nicks/${{repname}}.genrich.tn5nicks.bam\n" >> {output.tn5nicksFileList}
    echo -ne "${{repname}}\t{params.sample}\t{params.outdir}/${{repname}}.genrich.narrowPeak\n" >> {output.replicatePeakFileList}
done
echo -ne "{params.sample}\t{params.sample}\t{params.outdir}/{params.sample}.genrich.consensus.bed\n" > {output.consensusPeakFileList}

"""