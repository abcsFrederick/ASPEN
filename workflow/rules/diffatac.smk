# functions
def get_diffatac_input():
    expected_filelist = list()
    deffiles = list()
# diffatac only works for hg38 and mm10
    if GENOME == "mm10" or GENOME == "hg38":
        if CONTRASTS.shape[0] != 0:
            expected_filelist = [ join(RESULTSDIR,"peaks","genrich","DiffATAC","degs.done") ]
            expected_filelist.append([  join(RESULTSDIR,"peaks","genrich","DiffATAC","all_diff_atacs.html"),
                                    join(RESULTSDIR,"peaks","genrich","DiffATAC","all_diff_atacs.tsv")])
    return expected_filelist



rule atac_genrich_calculate_regions_of_interest_for_diffatac:
# """
# take genrich/macs consensus fixed width peaks and calculate ROIs
# """
    input:
        expand(join(RESULTSDIR,"peaks","genrich","fixed_width","{sample}.renormalized.fixed_width.consensus.narrowPeak"),sample=SAMPLES)
    output:
        roi              = join(RESULTSDIR,"peaks","genrich","fixed_width","ROI.narrowPeak"),
        roi_bed          = join(RESULTSDIR,"peaks","genrich","fixed_width","ROI.bed"),
        roi_gtf          = join(RESULTSDIR,"peaks","genrich","fixed_width","ROI.gtf"),
        roi_renormalized = join(RESULTSDIR,"peaks","genrich","fixed_width","ROI.renormalized.narrowPeak")
    params:
        script      =   "ccbr_atac_get_fixedwidth_peaks.sh",
        scriptsdir  =   SCRIPTSDIR,
        width       =   config["fixed_width"],
        roi_min_rep =   config['roi_min_replicates'],
        roi_min_spm =   config['roi_min_spm'],
    container: config["masterdocker"]
    shell:"""
set -exo pipefail
TMPDIR="/lscratch/$SLURM_JOB_ID"
if [ ! -d $TMPDIR ];then
    TMPDIR="/dev/shm"
fi

# get bedToGenePred genePredToGtf these are missing from the docker
cd $TMPDIR
wget -q -O bedToGenePred https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToGenePred
chmod a+x bedToGenePred
wget -q -O genePredToGtf https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToGtf
chmod a+x genePredToGtf

NPS=""
SAMPLENAMES=""
count=0
for i in {input}
do
    count=$((count+1))
    samplename=`basename $i | awk -F".renormalized" '{{print $1}}'`
    if [ "$count" == "1" ];then
        NPS="${{i}}"
        SAMPLENAMES="${{samplename}}"
        outdir=`dirname $i`
    else
        NPS="${{NPS}},${{i}}"
        SAMPLENAMES="${{SAMPLENAMES}},${{samplename}}"
    fi
done

cd $outdir
bash {params.scriptsdir}/{params.script} \
    --narrowpeaks $NPS \
    --replicatenames $SAMPLENAMES \
    --width {params.width} \
    --samplename "ROI" \
    --consensusnp {output.roi} \
    --consensusrenormnp {output.roi_renormalized} \
    --consensusminrep {params.roi_min_rep} \
    --consensusminspm {params.roi_min_spm} \
    --tmpdir $TMPDIR \
    --scriptsfolder {params.scriptsdir}
cut -f1-6 {output.roi} > {output.roi_bed}
bedSort {output.roi_bed} {output.roi_bed}
${{TMPDIR}}/bedToGenePred {output.roi_bed} /dev/stdout | ${{TMPDIR}}/genePredToGtf file /dev/stdin {output.roi_gtf}
"""

rule get_counts_table:
    input:
        files=expand(join(RESULTSDIR,"peaks","genrich","{sample}.genrich.tn5nicksbedfiles"),sample=SAMPLES),
        gtf=rules.atac_genrich_calculate_regions_of_interest_for_diffatac.output.roi_gtf
    output:
        counts=join(RESULTSDIR,"peaks","genrich","ROI.counts.tsv")
    params:
        scriptsdir  =   SCRIPTSDIR,
    threads: getthreads("get_counts_table")
    container: config['featurecountsdocker']
    shell:"""
set -exo pipefail
TMPDIR="/lscratch/$SLURM_JOB_ID"
if [ ! -d $TMPDIR ];then
    TMPDIR="/dev/shm"
fi
cd $TMPDIR
unset PYTHONPATH

bams=""
count=0
for f in {input.files}
do
while read a b c;do
    count=$((count+1))
    if [ "$count" == "1" ];then
        bams="$c"
    else
        bams="$bams $c"
    fi
done < $f
done

featureCounts -a {input.gtf} -o ${{TMPDIR}}/`basename {output.counts}` -s 0 -t transcript -g transcript_id -T {threads} $bams > featureCounts.log 2>&1
cat ${{TMPDIR}}/`basename {output.counts}` | python {params.scriptsdir}/_featureCounts_header_fix.py - > {output.counts}
"""

rule diffatac:
    input:
        counts      = rules.get_counts_table.output.counts,
    output:
        degsdone    = join(RESULTSDIR,"peaks","genrich","DiffATAC","degs.done"),
        sampleinfo  = join(RESULTSDIR,"peaks","genrich","DiffATAC","sampleinfo.txt"),
    params:
        contrasts   = config['contrasts'],
        scriptsdir  = SCRIPTSDIR,
        genome      = config['genome'],
        fc_cutoff   = config['contrasts_fc_cutoff'],
        fdr_cutoff  = config['contrasts_fdr_cutoff'],
        manifest    = config['samplemanifest'],
    container: config['baser']
    shell:"""
set -exo pipefail
TMPDIR="/lscratch/$SLURM_JOB_ID"
if [ ! -d $TMPDIR ];then
    TMPDIR="/dev/shm"
fi
outdir=$(dirname {output.degsdone})
if [ ! -d $outdir ];then mkdir $outdir;fi
cd $outdir
tail -n +2 {params.manifest} | cut -f1-2 | sort > {output.sampleinfo}

cat {params.contrasts}
while read g1 g2;do
    echo "g1=$g1"
    echo "g2=$g2"
    Rscript {params.scriptsdir}/DESeq2_runner.R \\
        --countsmatrix {input.counts} \\
        --genome {params.genome} \\
        --coldata {output.sampleinfo} \\
        --contrastnumerator $g1 \\
        --contrastdenominator $g2 \\
        --foldchange {params.fc_cutoff} \\
        --fdr {params.fdr_cutoff} \\
        --indexcols Geneid \\
        --excludecols Chr,Start,End,Strand,Length \\
        --outdir $outdir \\
        --tmpdir $TMPDIR \\
        --scriptsdir {params.scriptsdir}
done < {params.contrasts} && \
touch {output.degsdone}

"""

rule diffatac_aggregate:
    input:
        counts      = rules.get_counts_table.output.counts,
        degsdone    = rules.diffatac.output.degsdone,
        sampleinfo  = rules.diffatac.output.sampleinfo,
    output:
        alldegshtml = join(RESULTSDIR,"peaks","genrich","DiffATAC","all_diff_atacs.html"),
        alldegstsv  = join(RESULTSDIR,"peaks","genrich","DiffATAC","all_diff_atacs.tsv"),
    params:
        contrasts   = config['contrasts'],
        scriptsdir  = SCRIPTSDIR,
        genome      = config['genome'],
        fc_cutoff   = config['contrasts_fc_cutoff'],
        fdr_cutoff  = config['contrasts_fdr_cutoff'],
    container:config['baser']
    shell:"""
set -exo pipefail
TMPDIR="/lscratch/$SLURM_JOB_ID"
if [ ! -d $TMPDIR ];then
    TMPDIR="/dev/shm"
fi
outdir=$(dirname {input.degsdone})
cd $outdir
Rscript {params.scriptsdir}/aggregate_results_runner.R \\
    --countsmatrix {input.counts} \\
    --diffatacdir $outdir \\
    --coldata {input.sampleinfo} \\
    --foldchange {params.fc_cutoff} \\
    --fdr {params.fdr_cutoff} \\
    --indexcols Geneid \\
    --excludecols Chr,Start,End,Strand,Length \\
    --diffatacdir $outdir \\
    --tmpdir $TMPDIR \\
    --scriptsdir {params.scriptsdir}
"""
