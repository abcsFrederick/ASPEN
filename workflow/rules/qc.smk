rule fastqc:
# """
# Run FASTQC on:
# * Raw fastqs
# * Blacklist filtered trimmed fastqs
# """
    input:
        expand(join(WORKDIR,"fastqs","{replicate}.R1.fastq.gz"),replicate=REPLICATES),
        expand(join(WORKDIR,"fastqs","{replicate}.R2.fastq.gz"),replicate=REPLICATES),
        expand(rules.remove_BL.output.R1,replicate=REPLICATES),
        expand(rules.remove_BL.output.R2,replicate=REPLICATES),
    output:
        expand(join(QCDIR,"fastqc","{replicate}.R1_fastqc.zip"), replicate=REPLICATES),
        expand(join(QCDIR,"fastqc","{replicate}.R2_fastqc.zip"), replicate=REPLICATES),
        expand(join(QCDIR,"fastqc","{replicate}.R1.noBL_fastqc.zip"), replicate=REPLICATES),
        expand(join(QCDIR,"fastqc","{replicate}.R2.noBL_fastqc.zip"), replicate=REPLICATES),
    params:
        outdir=join(QCDIR,"fastqc"),
    threads: getthreads("fastqc")
    container: config["masterdocker"]
    shell: """
    fastqc {input} -t {threads} -o {params.outdir};
    """      

#########################################################

rule atac_tss:
# """
# * Creates read density profile near known TSS sites.
# * hg38/hg19/mm10 are the genomes supported, ie., TSS information 
# is deposited in the resources/db folder of the pipeline
# * per-replicate "tss.txt" file is generated --> read by MultiQC later
# """
    input:
        tagalign=join(RESULTSDIR,"tagAlign","{replicate}.tagAlign.gz")
    output:
        tss=join(QCDIR,"tss","{replicate}.tss.txt")
    params:
        replicate="{replicate}",
        workdir=RESULTSDIR,
        qcdir=QCDIR,
        indexdir=INDEXDIR,
        tssbed=TSSBED,
        genome=GENOME,
        scriptsdir=SCRIPTSDIR,
        script="ccbr_atac_tagAlign2TSS.bash",
    container: config["masterdocker"]
    threads: getthreads("atac_tss")    
    shell:"""
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then cd /lscratch/${{SLURM_JOB_ID}};else cd /dev/shm;fi
tagAlign=$(basename {input.tagalign})
outfn=$(basename {output.tss})

rsync -az --progress --verbose {input.tagalign} $tagAlign

bash {params.scriptsdir}/{params.script} \
--tagaligngz $tagAlign \
--genome {params.genome} \
--tsstxt $outfn \
--ncpus {threads} \
--tsstargz {params.tssbed}

rsync -az --progress --verbose $outfn {output.tss} && rm -f *
"""

#########################################################

rule atac_fld:
# """
# Calculate the Fragment-Length-Distribution per replicate
# Output: "fld.txt" file per-replicate --> used by MultiQC later
# """
    input:
        dedupbam=rules.align.output.dedupbam
    output:
        fld=join(QCDIR,"fld","{replicate}.fld.txt")
    params:
        replicate="{replicate}",
        scriptsdir=SCRIPTSDIR,
        script="ccbr_atac_bam2FLD.py",
    container: config["masterdocker"]    
    shell:"""
python {params.scriptsdir}/{params.script} -i {input.dedupbam} -o {output.fld}
"""

#########################################################

rule jaccard:
# """
# Use bedtools jaccard metric to quantify pair-wise score between
# * replicate peaks using the same peak caller (macs2 or genrich)
# * consensus peaks using the same peak caller (macs2 or genrich)
# * replicates and consensus peaks using the same peak caller (macs2 or genrich)
# * replicates only, consensus only and replicates+consensus peaks combined for both
#   peak callers together
# Output:
# * "pairwise.txt" files with jaccard score for all possible pairs
# * "jaccard.pca.html" --> calculate PCA using "pairwise.txt" file and plot using
#   plotly to a HTML document
# """
    input:
        expand(join(RESULTSDIR,"peaks","macs2","{sample}.consensus.macs2.peakfiles"),sample=SAMPLES),
        expand(join(RESULTSDIR,"peaks","macs2","{sample}.replicate.macs2.peakfiles"),sample=SAMPLES),
        expand(join(RESULTSDIR,"peaks","genrich","{sample}.consensus.genrich.peakfiles"),sample=SAMPLES),
        expand(join(RESULTSDIR,"peaks","genrich","{sample}.replicate.genrich.peakfiles"),sample=SAMPLES),
    output:
        macs2_per_replicate_jaccardpca=join(QCDIR,"jaccard","macs2.replicate.jaccard.pca.html"),
        macs2_per_sample_jaccardpca=join(QCDIR,"jaccard","macs2.consensus.jaccard.pca.html"),
        macs2_per_consensus_replicate_jaccardpca=join(QCDIR,"jaccard","macs2.consensus_replicate.jaccard.pca.html"),
        genrich_per_replicate_jaccardpca=join(QCDIR,"jaccard","genrich.replicate.jaccard.pca.html"),
        genrich_per_sample_jaccardpca=join(QCDIR,"jaccard","genrich.consensus.jaccard.pca.html"),
        genrich_per_consensus_replicate_jaccardpca=join(QCDIR,"jaccard","genrich.consensus_replicate.jaccard.pca.html"),
        allmethods_per_replicate_jaccardpca=join(QCDIR,"jaccard","allmethods.replicate.jaccard.pca.html"),
        allmethods_per_sample_jaccardpca=join(QCDIR,"jaccard","allmethods.consensus.jaccard.pca.html"),
        allmethods_per_consensus_replicate_jaccardpca=join(QCDIR,"jaccard","allmethods.consensus_replicate.jaccard.pca.html")
    params:
        workdir=RESULTSDIR,
        qcdir=QCDIR,
        indexdir=INDEXDIR,
        scriptsdir=SCRIPTSDIR,
        genome=GENOME,
        script="ccbr_jaccard_pca.bash",
    container: config["masterdocker"]
    threads: getthreads("jaccard")
    shell:"""
set -e -x -o pipefail
if [ -w "/lscratch/${{SLURM_JOB_ID}}" ];then TMPDIR="/lscratch/${{SLURM_JOB_ID}}";else cd TMPDIR="/dev/shm";fi
for f in {input};do
    rsync -Laz --progress $f $TMPDIR/
    while read replicateName sampleName file;do
        rsync -Laz --progress $file $TMPDIR/
    done < $f
done
cd $TMPDIR

for f in $(ls *narrowPeak);do
    bedSort $f $f
done

min_peaks=1000

awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *consensus.macs2.peakfiles > macs2.consensus.peakfiles.tmp
awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *replicate.macs2.peakfiles > macs2.replicate.peakfiles.tmp
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "$min_peaks" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < macs2.replicate.peakfiles.tmp > macs2.replicate.peakfiles
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "$min_peaks" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < macs2.consensus.peakfiles.tmp > macs2.consensus.peakfiles
cat macs2.consensus.peakfiles macs2.replicate.peakfiles > macs2.consensus_replicate.peakfiles
awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *consensus.genrich.peakfiles > genrich.consensus.peakfiles.tmp
awk -F"/" -v OFS="" "{{print \$1,\$NF}}" *replicate.genrich.peakfiles > genrich.replicate.peakfiles.tmp
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "$min_peaks" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < genrich.replicate.peakfiles.tmp > genrich.replicate.peakfiles
while read x y a;do w=$(wc -l $a|awk "{{print \$1}}"); if [ "$w" -gt "$min_peaks" ]; then echo -ne "$x\t$y\t$a\n" ;fi;done < genrich.consensus.peakfiles.tmp > genrich.consensus.peakfiles
cat genrich.consensus.peakfiles genrich.replicate.peakfiles > genrich.consensus_replicate.peakfiles

for m in "macs2" "genrich";do
    while read replicateName sampleName file;do
        echo -ne "${{replicateName}}_${{m}}\t${{sampleName}}_${{m}}\t${{file}}\n"
    done < ${{m}}.replicate.peakfiles
done > allmethods.replicate.peakfiles
for m in "macs2" "genrich";do
    while read replicateName sampleName file;do
        echo -ne "${{replicateName}}_${{m}}\t${{sampleName}}_${{m}}\t${{file}}\n"
    done < ${{m}}.consensus.peakfiles
done > allmethods.consensus.peakfiles
for m in "macs2" "genrich";do
    while read replicateName sampleName file;do
        echo -ne "${{replicateName}}_${{m}}\t${{sampleName}}_${{m}}\t${{file}}\n"
    done < ${{m}}.consensus_replicate.peakfiles
done > allmethods.consensus_replicate.peakfiles

for m in "macs2" "genrich" "allmethods";do
    for f in "consensus" "replicate" "consensus_replicate";do
        bash  {params.scriptsdir}/{params.script} \
        --inputfilelist ${{m}}.${{f}}.peakfiles \
        --pairwise ${{m}}.${{f}}.jaccard.pairwise.txt \
        --pcahtml ${{m}}.${{f}}.jaccard.pca.html \
        --scriptsfolder {params.scriptsdir}
        rsync -az --progress --verbose --remove-source-files ${{m}}.${{f}}.jaccard.pairwise.txt {params.qcdir}/jaccard/
        rsync -az --progress --verbose --remove-source-files ${{m}}.${{f}}.jaccard.pca.html {params.qcdir}/jaccard/
    done
done
"""

#########################################################