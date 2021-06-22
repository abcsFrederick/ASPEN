#!/bin/bash
set -e -x -o pipefail

callGenrichPeaks(){
BAM=$1
REPNAME=$2 
PEAKFILE=$3
GENOME=$4
GENOMEFILE=$5
EXCLUDELIST=$6

READSBEDFILE=${REPNAME}.reads.bed
READSBWFILE=${REPNAME}.bw
NICKSBEDFILE=${REPNAME}.tn5nicks.bed
NICKSBAMFILE=${REPNAME}.tn5nicks.bam 
#Genrich -t $BAM -o $PEAKFILE -j -r -e $EXCLUDELIST -v -s 5 -m 6 -b $READSBEDFILE -q 1 -l 200 -g 200
Genrich -t $BAM -o $PEAKFILE -j -r -e $EXCLUDELIST -v \
 -s $GENRICH_S \
 -m $GENRICH_M \
 -b $READSBEDFILE \
 -q $GENRICH_Q \
 -l $GENRICH_L \
 -g $GENRICH_G
processReadsBed $READSBEDFILE $READSBWFILE $GENOME $GENOMEFILE $NICKSBEDFILE $NICKSBAMFILE
}

processReadsBed() {
BED=$1
BW=$2
GENOME=$3
GENOMEFILE=$4
NICKSBED=$5
NICKSBAM=$6
if [ "$GENOME" == "hg19" ]; then
effectiveSize=2700000000
elif [ "$GENOME" == "hg38" ]; then
effectiveSize=2700000000
elif [ "$GENOME" == "mm9" ]; then
effectiveSize=2400000000
elif [ "$GENOME" == "mm10" ]; then
effectiveSize=2400000000
fi
bedSort $BED $BED
sf=$(awk -F"\t" -v size=$effectiveSize '{sum=sum+$3-$2}END{print sum/size}' $BED)
genomeCoverageBed -bg -i $BED -g $GENOMEFILE |awk -F"\t" -v OFS="\t" -v sf=$sf '{print $1,$2,$3,$4/sf}' > ${BED}.bg
bedGraphToBigWig ${BED}.bg $GENOMEFILE $BW
mv $BW ${OUTDIR}/bigwig/

awk -F"\t" -v OFS="\t" '{if ($3-$2==100) {print $1,$2+50,$2+51} }' $BED > $NICKSBED
awk -F"\t" -v OFS="\t" '{if ($3-$2!=100) {print $1,$2+50,$2+51} }' $BED >> $NICKSBED
awk -F"\t" -v OFS="\t" '{if ($3-$2!=100) {print $1,$3-51,$3-50} }' $BED >> $NICKSBED
bedSort $NICKSBED $NICKSBED
awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $NICKSBED > ${NICKSBED%.*}.tmp.bed
bedToBam -i ${NICKSBED%.*}.tmp.bed -g $GENOMEFILE > $NICKSBAM
samtools sort -@4 -o ${NICKSBAM%.*}.sorted.bam $NICKSBAM
mv ${NICKSBAM%.*}.sorted.bam $NICKSBAM
rm -f ${NICKSBED%.*}.tmp.bed ${BED}.bg
rm -f $BED
pigz -f -p4 $NICKSBED
mv ${NICKSBED}.gz ${OUTDIR}/tn5nicks/
mv $NICKSBAM ${OUTDIR}/tn5nicks/
samtools index ${OUTDIR}/tn5nicks/${NICKSBAM}
}

ARGPARSE_DESCRIPTION="call atac-seq peaks using genrich v0.6 with 2 replicates" 
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--bamfiles',required=True, nargs = '+', help = 'space separated list of BAM files')
parser.add_argument('--outdir',required=True, help= 'outputfolder')

parser.add_argument('--genome',required=True,help="hg19/38 or mm9/10")
parser.add_argument('--genomefile',required=True,help=".genome file")
parser.add_argument('--samplename',required=True, help='samplename,i.e., all replicates belong to this sample')

parser.add_argument('--genrich_s',required=False,default=5,help="Genrich -s parameter")
parser.add_argument('--genrich_m',required=False,default=6,help="Genrich -m parameter")
parser.add_argument('--genrich_q',required=False,default=1,help="Genrich -q parameter")
parser.add_argument('--genrich_l',required=False,default=200,help="Genrich -l parameter")
parser.add_argument('--genrich_g',required=False,default=200,help="Genrich -g parameter")


# parser.add_argument('--pooledpeakfile',required=False, help='output narrowPeak file for both replicates combined')
# parser.add_argument('--consensusbedfile',required=False, help='consensus overlapping MAX majority peaks in bed format')

# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="True", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=0.693147, help='default qfiltering value is 0.693147 (-log10 of 0.5) for q=0.5')

parser.add_argument('--scriptsfolder',required=False, default='/opt2', help='folder where the scripts are... used for debuging without rebuilding the docker')
EOF

cd $OUTDIR

CONSENSUSBEDFILE="${SAMPLENAME}.genrich.consensus.bed"
POOLEDPEAKFILE="${SAMPLENAME}.genrich.pooled.narrowPeak"

nreplicates=${#BAMFILES[@]}
BAMREP1=${BAMFILES[0]}
REP1NAME=`echo $(basename $BAMREP1)|awk -F".bam" '{print $1}'`
PEAKFILE1="${REP1NAME}.genrich.narrowPeak"
if [ "$nreplicates" -ge 2 ]; then
	BAMREP2=${BAMFILES[1]}
	REP2NAME=`echo $(basename $BAMREP2)|awk -F".bam" '{print $1}'`
	PEAKFILE2="${REP2NAME}.genrich.narrowPeak"
fi	
if [ "$nreplicates" -ge 3 ]; then
	BAMREP3=${BAMFILES[2]}
	REP3NAME=`echo $(basename $BAMREP3)|awk -F".bam" '{print $1}'`
	PEAKFILE2="${REP2NAME}.genrich.narrowPeak"
fi	
if [ "$nreplicates" -ge 4 ]; then
	BAMREP4=${BAMFILES[3]}
	REP4NAME=`echo $(basename $BAMREP4)|awk -F".bam" '{print $1}'`
	PEAKFILE4="${REP4NAME}.genrich.narrowPeak"
fi	


excludelist=$(samtools view -H $BAMREP1|grep ^@SQ|cut -f2|sed "s/SN://g"|awk "\$1 !~ /[0-9]$/"|grep -v "chrX"|grep -v "chrY"|tr "\n" ","|awk "{print substr(\$_,1,length(\$_)-1)}")
# echo $excludelist

# replicate 1 peak calling
callGenrichPeaks $BAMREP1 $REP1NAME $PEAKFILE1 $GENOME $GENOMEFILE "$excludelist"
# 
# replicate 2 peak calling
if [ "$nreplicates" -ge 2 ]; then
callGenrichPeaks $BAMREP2 $REP2NAME $PEAKFILE2 $GENOME $GENOMEFILE "$excludelist"
fi
# 
# replicate 3 peak calling
if [ "$nreplicates" -ge 3 ]; then
callGenrichPeaks $BAMREP3 $REP3NAME $PEAKFILE3 $GENOME $GENOMEFILE "$excludelist"
fi
# 
# replicate 4 peak calling
if [ "$nreplicates" -ge 4 ]; then
callGenrichPeaks $BAMREP4 $REP4NAME $PEAKFILE4 $GENOME $GENOMEFILE "$excludelist"
fi


if [ "$nreplicates" -ge 2 ]; then

	# merged peak calling and consensus peak calling
	READSBEDFILE=${POOLEDPEAKFILE%.*}.reads.bed
	READSBWFILE=${POOLEDPEAKFILE%.*}.reads.bw
	NICKSBEDFILE=${POOLEDPEAKFILE%.*}.tn5nicks.bed
	NICKSBAMFILE=${POOLEDPEAKFILE%.*}.tn5nicks.bam

	if [ "$nreplicates" -eq "2" ];then
	Genrich  -t $BAMREP1,$BAMREP2 -o $POOLEDPEAKFILE  -j -r -e $excludelist -v -s $GENRICH_S -m $GENRICH_M -b $READSBEDFILE -q $GENRICH_Q -l $GENRICH_L -g $GENRICH_G
	python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --qfilter $QFILTER --peakfiles $PEAKFILE1 $PEAKFILE2 --outbed $CONSENSUSBEDFILE
	elif [ "$nreplicates" -eq "3" ];then
	Genrich  -t $BAMREP1,$BAMREP2,$BAMREP3 -o $POOLEDPEAKFILE  -j -r -e $excludelist -v -s $GENRICH_S -m $GENRICH_M -b $READSBEDFILE -q $GENRICH_Q -l $GENRICH_L -g $GENRICH_G
	python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --qfilter $QFILTER --peakfiles $PEAKFILE1 $PEAKFILE2 $PEAKFILE3 --outbed $CONSENSUSBEDFILE
	elif [ "$nreplicates" -eq "4" ];then
	Genrich  -t $BAMREP1,$BAMREP2,$BAMREP3,$BAMREP4 -o $POOLEDPEAKFILE  -j -r -e $excludelist -v -s $GENRICH_S -m $GENRICH_M -b $READSBEDFILE -q $GENRICH_Q -l $GENRICH_L -g $GENRICH_G
	python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --qfilter $QFILTER --peakfiles $PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $PEAKFILE4 --outbed $CONSENSUSBEDFILE
	fi


	processReadsBed $READSBEDFILE $READSBWFILE $GENOME $GENOMEFILE $NICKSBEDFILE $NICKSBAMFILE
fi

if [ "$nreplicates" -eq 1 ];then
	for f in `ls ${REP1NAME}.genrich*`;do
		newf=$(echo $f|sed "s/.genrich./.genrich.pooled./g")
		cp $f $newf
	done
	cut -f1-3 $PEAKFILE1 > $CONSENSUSBEDFILE
fi

# if [ $CONSENSUSBEDFILE ];then
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


files="$PEAKFILE1"
if [ "$nreplicates" -eq "2" ];then files="$PEAKFILE1 $PEAKFILE2 $POOLEDPEAKFILE"; fi
if [ "$nreplicates" -eq "3" ];then files="$PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $POOLEDPEAKFILE"; fi
if [ "$nreplicates" -eq "4" ];then files="$PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $PEAKFILE4 $POOLEDPEAKFILE"; fi


if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  for f in $files;do
	filteredpeakfile=$(echo $f|sed "s/.narrowPeak/.qfilter.narrowPeak/g")
	awk -F"\t" -v q=$qvalue "{if (\$9>q){print}}" $f > $filteredpeakfile
	npeaks=$(wc -l $f|awk '{print $1}')
	if [ "$npeaks" -gt "0" ];then
		Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $f -a ${f}.annotated -g $GENOME -l ${f}.genelist -f ${f}.annotation_summary
		cut -f1,2 ${f}.annotation_summary > ${f}.annotation_distribution
	else
		touch ${f}.annotated
		touch ${f}.genelist
		touch ${f}.annotation_summary
		touch ${f}.annotation_distribution
	fi
	npeaks_filtered=$(wc -l $filteredpeakfile|awk '{print $1}')
	if [ "$npeaks_filtered" -gt "0" ];then
		Rscript ${SCRIPTSFOLDER}/ccbr_annotate_peaks.R -n $filteredpeakfile -a ${filteredpeakfile}.annotated -g $GENOME -l ${filteredpeakfile}.genelist -f ${filteredpeakfile}.annotation_summary
		cut -f1,2 ${filteredpeakfile}.annotation_summary > ${filteredpeakfile}.annotation_distribution
	else
		touch ${filteredpeakfile}.annotated
		touch ${filteredpeakfile}.genelist
		touch ${filteredpeakfile}.annotation_summary
		touch ${filteredpeakfile}.annotation_distribution
	fi		
  done
fi
