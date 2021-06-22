def get_sample_tagAlignFiles_for_group(wildcards):
    samples=GROUP2SAMPLES[wildcards.grp]
    d=dict()
    for i,s in enumerate(samples):
        d["tagAlign"+str(i+1)]=join(RESULTSDIR,"tagAlign",s+".tagAlign.gz")
    return d

rule atac_macs_peakcalling:
    input:
        unpack(get_sample_tagAlignFiles_for_group)
    params:
        genome=GENOME,
        genomefile=GENOMEFILE,
        outdir=join(RESULTSDIR,"peaks","macs2"),
        scriptsdir=SCRIPTSDIR,
        qcdir=QCDIR,
        script="ccbr_atac_macs2_peak_calling.bash",
        grp="{grp}"
    output:
        consensusPeakFileList=join(RESULTSDIR,"peaks","macs2","{grp}.consensus.macs2.peakfiles"),
        replicatePeakFileList=join(RESULTSDIR,"peaks","macs2","{grp}.replicate.macs2.peakfiles"),
        tn5nicksFileList=join(RESULTSDIR,"peaks","macs2","{grp}.macs2.tn5nicksbedfiles")
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
    --samplename {params.grp} \
    --scriptsfolder {params.scriptsdir} \
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
    repname=`basename $f|awk -F".tagAlign" '{{print $1}}'`
    echo -ne "${{repname}}\t{params.grp}\t{params.outdir}/tn5nicks/${{repname}}.macs2.tn5nicks.bam\n" >> {output.tn5nicksFileList}
    echo -ne "${{repname}}\t{params.grp}\t{params.outdir}/${{repname}}.macs2.narrowPeak\n" >> {output.replicatePeakFileList}
done
echo -ne "{params.grp}\t{params.grp}\t{params.outdir}/{params.grp}.macs2.consensus.bed\n" > {output.consensusPeakFileList}

"""