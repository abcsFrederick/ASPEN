#!/bin/bash
set -e -x -o pipefail

callGenrichPeaks(){
BAM=$1
PEAKFILE=$2
EXCLUDELIST=$3
READSBEDFILE=$4

Genrich -t $BAM \
 -o $PEAKFILE \
 -j \
 -r \
 -v \
 -e $EXCLUDELIST \
 -s $GENRICH_S \
 -m $GENRICH_M \
 -b $READSBEDFILE \
 -q $GENRICH_Q \
 -l $GENRICH_L \
 -g $GENRICH_G \
 -d $GENRICH_D
}

callGenrichPeaksProcessReadsBed(){
BAM=$1
REPNAME=$2 
PEAKFILE=$3
EXCLUDELIST=$4

READSBEDFILE=${REPNAME}.reads.bed
READSBWFILE=${REPNAME}.bw
NICKSBEDFILE=${REPNAME}.tn5nicks.bed
NICKSBAMFILE=${REPNAME}.tn5nicks.bam 
#Genrich -t $BAM -o $PEAKFILE -j -r -e $EXCLUDELIST -v -s 5 -m 6 -b $READSBEDFILE -q 1 -l 200 -g 200
callGenrichPeaks $BAM $PEAKFILE $EXCLUDELIST $READSBEDFILE 
processReadsBed $READSBEDFILE $READSBWFILE $NICKSBEDFILE $NICKSBAMFILE
}

processReadsBed() {
BED=$1
BW=$2
NICKSBED=$3
NICKSBAM=$4
effectivegenomesize=$EFFECTIVEGENOMESIZE

bedSort $BED $BED
sf=$(awk -F"\t" -v size=$effectiveSize '{sum=sum+$3-$2}END{print sum/size}' $BED)
genomeCoverageBed -bg -i $BED -g $GENOMEFILE |awk -F"\t" -v OFS="\t" -v sf=$sf '{print $1,$2,$3,$4/sf}' > ${BED}.bg
# some coordinates go out of chromosome size in macs2 bam2bw conversion.... using a working around there which is copied here as well
bn=$(basename $BED)
awk -F"\t" -v OFS="\t" '{print $1,"0",$2}' $GENOMEFILE > /dev/shm/${bn}.genome.bed
bedtools intersect -a ${BED}.bg -b /dev/shm/${bn}.genome.bed > /dev/shm/${bn}.bg
mv /dev/shm/${bn}.bg ${BED}.bg
rm -f /dev/shm/${bn}.genome.bed
#
bedGraphToBigWig ${BED}.bg $GENOMEFILE $BW
mv $BW ${OUTDIR}/bigwig/
d=$GENRICH_D
e=$((d/2))
awk -F"\t" -v OFS="\t" -v d=d -v e=e '{if ($3-$2==d) {print $1,$2+e,$2+e+1} }' $BED > $NICKSBED
awk -F"\t" -v OFS="\t" -v d=d -v e=e '{if ($3-$2!=d) {print $1,$2+e,$2+e+1} }' $BED >> $NICKSBED
awk -F"\t" -v OFS="\t" -v d=d -v e=e '{if ($3-$2!=d) {print $1,$3-e-1,$3-e} }' $BED >> $NICKSBED
bedSort $NICKSBED $NICKSBED
awk -F"\t" -v OFS="\t" '{print $1,$2,$3,$1":"$2"-"$3}' $NICKSBED > ${NICKSBED%.*}.tmp.bed
bedToBam -i ${NICKSBED%.*}.tmp.bed -g $GENOMEFILE > $NICKSBAM
samtools sort -@4 -o ${NICKSBAM%.*}.sorted.bam $NICKSBAM
mv ${NICKSBAM%.*}.sorted.bam $NICKSBAM
rm -f ${NICKSBED%.*}.tmp.bed ${BED}.bg
pigz -f -p4 $BED
mv ${BED}.gz $READSBEDFOLDER/
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
parser.add_argument('--effectivegenomesize',required=True,help="effective genome size")
parser.add_argument('--samplename',required=True, help='samplename,i.e., all replicates belong to this sample')

parser.add_argument('--genrich_d',required=False,default=100,help="Genrich -d parameter")
parser.add_argument('--genrich_s',required=False,default=5,help="Genrich -s parameter")
parser.add_argument('--genrich_m',required=False,default=6,help="Genrich -m parameter")
parser.add_argument('--genrich_q',required=False,default=1,help="Genrich -q parameter")
parser.add_argument('--genrich_l',required=False,default=100,help="Genrich -l parameter")
parser.add_argument('--genrich_g',required=False,default=100,help="Genrich -g parameter")
parser.add_argument('--runchipseeker',required=False,default="False",help="run chipseeker")
parser.add_argument('--readsbedfolder',required=True,default="False",help="folder to save reads used by Genrich")



# parser.add_argument('--pooledpeakfile',required=False, help='output narrowPeak file for both replicates combined')
# parser.add_argument('--consensusbedfile',required=False, help='consensus overlapping MAX majority peaks in bed format')

# only required if filtering
parser.add_argument('--filterpeaks',required=False, default="True", help='filterpeaks by qvalue: True or False')
parser.add_argument('--qfilter',required=False, default=1, help='default qfiltering value is 1 (-log10 of 0.1) for q=0.1')

parser.add_argument('--scriptsfolder',required=True, help='folder where the scripts are ... usually workflow/scripts')
EOF

cd $OUTDIR

CONSENSUSBEDFILE="${SAMPLENAME}.genrich.consensus.bed"
CONSENSUSFILTER=0.5 # more than 50% of the replicates must be represented in each peak.
POOLEDPEAKFILE="${SAMPLENAME}.genrich.pooled.narrowPeak"

nreplicates=${#BAMFILES[@]}
BAMREP1=${BAMFILES[0]}
REP1NAME=`echo $(basename $BAMREP1)|awk -F".qsorted.bam" '{print $1}'`
PEAKFILE1="${REP1NAME}.genrich.narrowPeak"
if [ "$nreplicates" -ge 2 ]; then
	BAMREP2=${BAMFILES[1]}
	REP2NAME=`echo $(basename $BAMREP2)|awk -F".qsorted.bam" '{print $1}'`
	PEAKFILE2="${REP2NAME}.genrich.narrowPeak"
fi	
if [ "$nreplicates" -ge 3 ]; then
	BAMREP3=${BAMFILES[2]}
	REP3NAME=`echo $(basename $BAMREP3)|awk -F".qsorted.bam" '{print $1}'`
	PEAKFILE3="${REP3NAME}.genrich.narrowPeak"
fi	
if [ "$nreplicates" -ge 4 ]; then
	BAMREP4=${BAMFILES[3]}
	REP4NAME=`echo $(basename $BAMREP4)|awk -F".qsorted.bam" '{print $1}'`
	PEAKFILE4="${REP4NAME}.genrich.narrowPeak"
fi	


excludelist=$(samtools view -H $BAMREP1|grep ^@SQ|cut -f2|sed "s/SN://g"|awk "\$1 !~ /[0-9]$/"|grep -v "chrX"|grep -v "chrY"|tr "\n" ","|awk "{print substr(\$_,1,length(\$_)-1)}")
# echo $excludelist

# replicate 1 peak calling
callGenrichPeaksProcessReadsBed $BAMREP1 $REP1NAME $PEAKFILE1 "$excludelist"
# 
# replicate 2 peak calling
if [ "$nreplicates" -ge 2 ]; then
callGenrichPeaksProcessReadsBed $BAMREP2 $REP2NAME $PEAKFILE2 "$excludelist"
fi
# 
# replicate 3 peak calling
if [ "$nreplicates" -ge 3 ]; then
callGenrichPeaksProcessReadsBed $BAMREP3 $REP3NAME $PEAKFILE3 "$excludelist"
fi
# 
# replicate 4 peak calling
if [ "$nreplicates" -ge 4 ]; then
callGenrichPeaksProcessReadsBed $BAMREP4 $REP4NAME $PEAKFILE4 "$excludelist"
fi


if [ "$nreplicates" -ge 2 ]; then

	# merged peak calling and consensus peak calling
	READSBEDFILE=${POOLEDPEAKFILE%.*}.reads.bed
	READSBWFILE=${POOLEDPEAKFILE%.*}.reads.bw
	NICKSBEDFILE=${POOLEDPEAKFILE%.*}.tn5nicks.bed
	NICKSBAMFILE=${POOLEDPEAKFILE%.*}.tn5nicks.bam

	if [ "$nreplicates" -eq "2" ];then
	callGenrichPeaks "$BAMREP1,$BAMREP2" $POOLEDPEAKFILE "$excludelist" $READSBEDFILE
	python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --filter $CONSENSUSFILTER --peakfiles $PEAKFILE1 $PEAKFILE2 --outbed $CONSENSUSBEDFILE
	elif [ "$nreplicates" -eq "3" ];then
	callGenrichPeaks  "$BAMREP1,$BAMREP2,$BAMREP3" $POOLEDPEAKFILE  "$excludelist" $READSBEDFILE
	python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --filter $CONSENSUSFILTER --peakfiles $PEAKFILE1 $PEAKFILE2 $PEAKFILE3 --outbed $CONSENSUSBEDFILE
	elif [ "$nreplicates" -eq "4" ];then
	callGenrichPeaks  "$BAMREP1,$BAMREP2,$BAMREP3,$BAMREP4" $POOLEDPEAKFILE  "$excludelist" $READSBEDFILE
	python ${SCRIPTSFOLDER}/ccbr_get_consensus_peaks.py --filter $CONSENSUSFILTER --peakfiles $PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $PEAKFILE4 --outbed $CONSENSUSBEDFILE
	fi

	processReadsBed $READSBEDFILE $READSBWFILE $NICKSBEDFILE $NICKSBAMFILE
fi

if [ "$nreplicates" -eq 1 ];then
	for f in `ls ${REP1NAME}.genrich*`;do
		newf=$(echo $f|sed "s/.genrich./.genrich.pooled./g")
		cp $f $newf
	done
	cut -f1-3 $PEAKFILE1 > $CONSENSUSBEDFILE
fi

if [ "$RUNCHIPSEEKER" == "True" ];then
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
fi


files="$PEAKFILE1"
if [ "$nreplicates" -eq "2" ];then files="$PEAKFILE1 $PEAKFILE2 $POOLEDPEAKFILE"; fi
if [ "$nreplicates" -eq "3" ];then files="$PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $POOLEDPEAKFILE"; fi
if [ "$nreplicates" -eq "4" ];then files="$PEAKFILE1 $PEAKFILE2 $PEAKFILE3 $PEAKFILE4 $POOLEDPEAKFILE"; fi

#ChIPseeker in the container only works for hg19/hg38/mm10... so you cannot annotate other genomes here
genome_is_known=0
if [ "$RUNCHIPSEEKER" == "True" ];then
	if [ "$GENOME" == "hg19" ];then
		genome_is_known=1
	elif [ "$GENOME" == "hg38" ];then
		genome_is_known=1
	elif [ "$GENOME" == "mm10" ];then
		genome_is_known=1
	fi
	if [ "$genome_is_known" == "0" ];then
	RUNCHIPSEEKER="False"
	fi
fi

if [ $FILTERPEAKS == "True" ];then
  qvalue=$QFILTER
  for f in $files;do
	filteredpeakfile=$(echo $f|sed "s/.narrowPeak/.qfilter.narrowPeak/g")
	awk -F"\t" -v q=$qvalue "{if (\$9>q){print}}" $f > $filteredpeakfile
	if [ "$RUNCHIPSEEKER" == "True" ];then
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
	fi
  done
fi
