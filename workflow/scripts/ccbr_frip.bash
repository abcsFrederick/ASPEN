#!/bin/bash

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="fraction of reads in peaks (and DHS/Enhancers/Promoters [if genome is provided])" 
source /opt2/argparse.bash || \
source <(curl -s https://raw.githubusercontent.com/CCBR/Tools/master/scripts/argparse.bash) || \
exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--narrowpeak',required=True, help='input narrowPeak or bed file')
parser.add_argument('--peakcaller',required=False, default="", help='peakcaller')
parser.add_argument('--tagalign',required=True, help='input tagAlign.gz or bed.gz file')
parser.add_argument('--samplename',required=True, help='sample name')
parser.add_argument('--dhsbed',required=False, help='known DHS sites')
parser.add_argument('--promoterbed',required=False, help='known Promoter sites')
parser.add_argument('--enhancerbed',required=False, help='known Enhancer sites')
EOF

set -e -x -o pipefail
nreads=$(zcat $TAGALIGN|wc -l)

frip_reads=$(bedtools intersect -wa -a $TAGALIGN -b $NARROWPEAK|wc -l)
frip=$(echo "scale=3;$frip_reads/$nreads"|bc)
echo -ne "$SAMPLENAME\tFRiP${PEAKCALLER}\t$frip\n"

if [ $DHSBED ];then
if [ "$DHSBED" != "" ];then
fridhs_reads=$(bedtools intersect -wa -a $TAGALIGN -b $DHSBED|wc -l)
fridhs=$(echo "scale=3;$fridhs_reads/$nreads"|bc)
echo -ne "$SAMPLENAME\tFRiDHS\t$fridhs\n"
fi
fi

if [ $PROMOTERBED ];then
if [ "$PROMOTERBED" != "" ];then
fripro_reads=$(bedtools intersect -wa -a $TAGALIGN -b $PROMOTERBED|wc -l)
fripro=$(echo "scale=3;$fripro_reads/$nreads"|bc)
echo -ne "$SAMPLENAME\tFRiPro\t$fripro\n"
fi
fi

if [ $ENHANCERBED ];then
if [ "$ENHANCEDBED" != "" ];then
frienh_reads=$(bedtools intersect -wa -a $TAGALIGN -b $ENHANCERBED|wc -l)
frienh=$(echo "scale=3;$frienh_reads/$nreads"|bc)
echo -ne "$SAMPLENAME\tFRiEnh\t$frienh\n"
fi
fi
