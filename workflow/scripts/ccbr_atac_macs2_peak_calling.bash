#!/bin/bash

callPeaks(){
# inputs
	TAGALIGN=$1
	PREFIX=$2
	GENOME=$3
	SHIFTSIZE=$4
	EXTSIZE=$5
	FILTERPEAKS=$6
	QFILTER=$7
	GENOMEFILE=$8
	SCRIPTSFOLDER=$9
	OUTDIR=${10}

# call peaks
	if [ $GENOME == "hg19" ]; then g="hs";fi
	if [ $GENOME == "hg38" ]; then g="hs";fi
	if [ $GENOME == "mm10" ]; then g="mm";fi
	if [ $GENOME == "mm9" ]; then g="mm";fi
	macs2 callpeak -t $TAGALIGN -f BED -n $PREFIX -g $g -p 0.01 --shift -$SHIFTSIZE --extsize $EXTSIZE --keep-dup all --call-summits --nomodel --SPMR -B --outdir $OUTDIR

# remove duplicate peak calls
	mv ${PREFIX}_peaks.narrowPeak ${PREFIX}_peaks.narrowPeak.tmp
	sort -k9,9gr ${PREFIX}_peaks.narrowPeak.tmp|awk -F"\t" '!NF || !seen[$1":"$2"-"$3]++'|sort -k1,1 -k2,2n > ${PREFIX}_peaks.narrowPeak
	rm -f ${PREFIX}_peaks.narrowPeak.tmp
	mv ${PREFIX}_peaks.narrowPeak ${PREFIX}.narrowPeak

# qvalue filter of 0.5 ... very lenient
	if [ $FILTERPEAKS == "True" ];then
	  qvalue=$QFILTER
	  awk -F"\t" -v q=$qvalue '{if ($9>q){print}}' ${PREFIX}.narrowPeak > ${PREFIX}.qfilter.narrowPeak
	fi

# annotate
	for f in ${PREFIX}.narrowPeak ${PREFIX}.qfilter.narrowPeak;do
		npeaks=$(wc -l $f|awk '{print $1}')
		if [ "$npeaks" -gt "0" ];then
			Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $f  -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary 
			cut -f1,2 ${f}.annotation_summary > ${f}.annotation_distribution
		else
			touch ${f}.annotated
			touch ${f}.genelist
			touch ${f}.annotation_summary
			touch ${f}.annotation_distribution
		fi
	done

# save bigwig
# MAY NEED TO BE NORMALIZED FOR GENOME SIZE ... DONT KNOW FOR SURE
# MACS says
#
# MACS will save fragment pileup signal per million reads
#
# Genrich normalization is as per 1x genome coverage... and macs2 normalization is per million reads 
	bedSort ${PREFIX}_treat_pileup.bdg ${PREFIX}_treat_pileup.bdg
	bedGraphToBigWig ${PREFIX}_treat_pileup.bdg $GENOMEFILE ${PREFIX}.bw
	rm -f ${PREFIX}_treat_pileup.bdg
	rm -f ${PREFIX}_control_lambda.bdg
	mv ${PREFIX}.bw ${OUTDIR}/bigwig/${PREFIX}.bw

# save nicks bam and bed
  nicksBED=${PREFIX}.tn5nicks.bed
  nicksBAM=${PREFIX}.tn5nicks.bam
  zcat $TAGALIGN|awk -F"\t" -v OFS="\t" '{if ($6=="+") {print $1,$2,$2+1} else {print $1,$3,$3+1}}'> $nicksBED
  bedSort $nicksBED $nicksBED
  awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $nicksBED > ${nicksBED%.*}.tmp.bed
  bedToBam -i ${nicksBED%.*}.tmp.bed -g $GENOMEFILE > $nicksBAM
  samtools sort -@4 -o ${nicksBAM%.*}.sorted.bam $nicksBAM
  mv ${nicksBAM%.*}.sorted.bam $nicksBAM
  rm -f ${nicksBED%.*}.tmp.bed 
  pigz -f -p4 $nicksBED
  mv ${nicksBED}.gz ${OUTDIR}/tn5nicks/${nicksBED}.gz
  mv ${nicksBAM} ${OUTDIR}/tn5nicks/${nicksBAM}
  samtools index ${OUTDIR}/tn5nicks/${nicksBAM}
}

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="""This scripts takes in Tn5-shifted tagAligne files as input and
1. Calls atac-seq peaks using macs2
2. Filters low-quality script using q-value
3. Annotates peaks using ChIPseeker
4. Calls atac-seq peaks for pooled replicates and annotates them
5. Calls consensus peaks for multiple replicates and annotates them
* Currently supports upto 4 replicates
Output files:
* Per replicate:
	** narrowPeak file
	** annotated file from ChIPseeker
	** annotations summary file
	** annotation distribution file (for pie-charting)
	** genelist file
	** sorted tn5nicks BED file (in tn5nicks subfolder)
	** tn5nicks BAM file (required for downstream DiffBind analysis ... in tn5nicks subfolder)
	** bigwig files (in bigwig subfolder)
* Pooled narrowPeak file with annotations
* Consensus BED file with annotations 
"""
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--tagalignfiles',required=True, nargs = '+', help = 'space separated list of tagAlignGz files')
parser.add_argument('--outdir',required=True, help= 'outputfolder')

parser.add_argument('--genome',required=True,help="hg19/38 or mm9/10")
parser.add_argument('--genomefile',required=True,help=".genome file")

parser.add_argument('--extsize',required=False, default=200, help='extsize')
parser.add_argument('--shiftsize',required=False, default=100, help='shiftsize')
parser.add_argument('--samplename',required=True, help='samplename,i.e., all replicates belong to this sample')

# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="True", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=0.693147, help='default qfiltering value is 0.693147 (-log10 of 0.5) for q=0.5')

parser.add_argument('--scriptsfolder',required=True,  help='folder where the scripts are... probably <path to workflow>/scripts')

EOF

cd $OUTDIR

CONSENSUSBEDFILE="${SAMPLENAME}.macs2.consensus.bed"

nreplicates=${#TAGALIGNFILES[@]}
TAGALIGN1=${TAGALIGNFILES[0]}
REP1NAME=`echo $(basename $TAGALIGN1)|awk -F".tag" '{print $1}'`
if [ "$nreplicates" -ge 2 ]; then
	TAGALIGN2=${TAGALIGNFILES[1]}
	REP2NAME=`echo $(basename $TAGALIGN2)|awk -F".tag" '{print $1}'`
fi	
if [ "$nreplicates" -ge 3 ]; then
	TAGALIGN3=${TAGALIGNFILES[2]}
	REP3NAME=`echo $(basename $TAGALIGN3)|awk -F".tag" '{print $1}'`
fi	
if [ "$nreplicates" -ge 4 ]; then
	TAGALIGN4=${TAGALIGNFILES[3]}
	REP4NAME=`echo $(basename $TAGALIGN4)|awk -F".tag" '{print $1}'`
fi	


# replicate 1 peak calling
callPeaks $TAGALIGN1 ${REP1NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER $OUTDIR
# 
# replicate 2 peak calling
if [ "$nreplicates" -ge 2 ]; then
callPeaks $TAGALIGN2 ${REP2NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER $OUTDIR
fi
# 
# replicate 3 peak calling
if [ "$nreplicates" -ge 3 ]; then
callPeaks $TAGALIGN3 ${REP3NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER $OUTDIR
fi
# 
# replicate 4 peak calling
if [ "$nreplicates" -ge 4 ]; then
callPeaks $TAGALIGN4 ${REP4NAME}.macs2 $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER $OUTDIR
fi

if [ "$nreplicates" -eq 1 ];then
	for f in `ls ${REP1NAME}.macs2*`;do
		newf=$(echo $f|sed "s/.macs2./.macs2.pooled./g")
		cp $f $newf
	done
	cut -f1-3 ${REP1NAME}.macs2.narrowPeak > $CONSENSUSBEDFILE
fi

if [ "$nreplicates" -eq 2 ];then
pooled="${REP1NAME}_${REP2NAME}.pooled.tagAlign"
zcat $TAGALIGN1 $TAGALIGN2 > ${pooled}
fi

if [ "$nreplicates" -eq 3 ];then
pooled="${REP1NAME}_${REP2NAME}_${REP3NAME}.pooled.tagAlign"
zcat $TAGALIGN1 $TAGALIGN2 $TAGALIGN3 > ${pooled}
fi

if [ "$nreplicates" -eq 4 ];then
pooled="${REP1NAME}_${REP2NAME}_${REP3NAME}_${REP4NAME}.pooled.tagAlign"
zcat $TAGALIGN1 $TAGALIGN2 $TAGALIGN3 $TAGALIGN4 > ${pooled}
fi

# pooled peak calling
if [ "$nreplicates" -ge 2 ];then
bedSort ${pooled} ${pooled}
pigz -f -p4 ${pooled}
callPeaks ${pooled}.gz ${SAMPLENAME}.macs2.pooled $GENOME $SHIFTSIZE $EXTSIZE $FILTERPEAKS $QFILTER $GENOMEFILE $SCRIPTSFOLDER $OUTDIR
fi

# consensus peak calling
if [ "$nreplicates" -eq 2 ];then
python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --qfilter $QFILTER --peakfiles ${REP1NAME}.macs2.narrowPeak ${REP2NAME}.macs2.narrowPeak --outbed $CONSENSUSBEDFILE
fi

if [ "$nreplicates" -eq 3 ];then
python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --qfilter $QFILTER --peakfiles ${REP1NAME}.macs2.narrowPeak ${REP2NAME}.macs2.narrowPeak ${REP3NAME}.macs2.narrowPeak --outbed $CONSENSUSBEDFILE
fi

if [ "$nreplicates" -eq 4 ];then
python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --qfilter $QFILTER --peakfiles ${REP1NAME}.macs2.narrowPeak ${REP2NAME}.macs2.narrowPeak ${REP3NAME}.macs2.narrowPeak ${REP4NAME}.macs2.narrowPeak --outbed $CONSENSUSBEDFILE
fi

# if [ "$nreplicates" -ge 2 ];then
	f="$CONSENSUSBEDFILE"
	npeaks_consensus=$(wc -l $f|awk '{print $1}')
	if [ "$npeaks_consensus" -gt "0" ];then
		Rscript ${SCRIPTSFOLDER}/ccbr_annotate_bed.R -b $f -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary
		cut -f1,2 ${f}.annotation_summary > ${f}.annotation_distribution
	else
		touch ${f}.annotated
		touch ${f}.genelist
		touch ${f}.annotation_summary
		touch ${f}.annotation_distribution
	fi
# fi
