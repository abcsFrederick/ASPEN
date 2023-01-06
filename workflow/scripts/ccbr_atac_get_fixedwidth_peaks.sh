#!/usr/bin/env bash
set -e -x -o pipefail

ARGPARSE_DESCRIPTION="call atac-seq peaks using genrich v0.6 with 2 replicates" 
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--narrowpeaks',required=True, nargs = '+', help = 'space separated list of narrowPeak files')
parser.add_argument('--replicatenames',required=True, nargs = '+', help = 'space separated list of replicate names')
parser.add_argument('--width',required=False,default=500,help="default 500")
parser.add_argument('--samplename',required=True, help='samplename,i.e., all replicates belong to this sample')

parser.add_argument('--consensusnp',required=True,help="consensus narrowPeak file")
parser.add_argument('--consensusrenormnp',required=True,help="consensus renormalized narrowPeak file")

parser.add_argument('--tmpdir',required=False,default="/dev/shm",help="temp dir")

parser.add_argument('--scriptsfolder',required=True, help='folder where the scripts are ... usually workflow/scripts')
EOF

# nreplicates=${#NARROWPEAKS[@]}

replicateNames=$(echo $REPLICATENAMES|sed "s/ /,/g")

count=0
for np in NARROWPEAKS {
    count=$((count+1))
    outnp=$(echo $np|sed "s/.narrowPeak/.fixed_width.narrowPeak/g")
    if [ "$count" == "1" ];then
        OUTNPS="$outnp"
        # REPNAMES=$(basename $np | awk -F".narrowPeak" '{print $1}')
    else
        OUTNPS="$OUTNPS,$outnp"
        # repname=$(basename $np | awk -F".narrowPeak" '{print $1}')
        # REPNAMES="$REPNAMES,$repname"
    fi
    Rscript ${SCRIPTSFOLDER}/narrowPeak_to_fixed_width_peakSet.R \
    -i $np \
    -o $outnp \
    -w $WIDTH \
    -t $TMPDIR
}

Rscript ${SCRIPTSFOLDER}/fixed_width_peakSets_to_consensus_peakSet.R \
    -i $OUTNPS \
    -p $replicateNames \
    -o $CONSENSUSNP \
    -t $TMPDIR \
    -m 2 \
    -c 5

Rscript ${SCRIPTSFOLDER}/narrowPeak_normalize_pvalues.R \
    -i $CONSENSUSNP \
    -o $CONSENSUSRENORMNP