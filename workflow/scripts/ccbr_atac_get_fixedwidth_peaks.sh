#!/usr/bin/env bash
set -e -x -o pipefail

ARGPARSE_DESCRIPTION="call atac-seq peaks using genrich v0.6 with 2 replicates" 
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--narrowpeaks',required=True, nargs = '+', help = 'comma separated list of narrowPeak files')
parser.add_argument('--replicatenames',required=True, nargs = '+', help = 'comma separated list of replicate names')
parser.add_argument('--width',required=False,default=500,help="default 500")
parser.add_argument('--samplename',required=True, help='samplename,i.e., all replicates belong to this sample')

parser.add_argument('--consensusnp',required=True,help="consensus narrowPeak file")
parser.add_argument('--consensusrenormnp',required=True,help="consensus renormalized narrowPeak file")
parser.add_argument('--consensusminrep',required=False,default=2,type=int,help="minimum replicates that need to overlap for any peak to be considered a consensus peak")
parser.add_argument('--consensusminspm',required=False,default=5,type=int,help="cutoff score per million or normalized p-value. All minimum number of replicate peaks that need to overlap with a consensus peak should also be greater than this cutoff score per million or SPM")

parser.add_argument('--tmpdir',required=False,default="/dev/shm",help="temp dir")

parser.add_argument('--scriptsfolder',required=True, help='folder where the scripts are ... usually workflow/scripts')
EOF

# nreplicates=${#NARROWPEAKS[@]}

# REPLICATENAMES=$(echo $REPLICATENAMES|sed "s/,/ /g")
NARROWPEAKS=$(echo $NARROWPEAKS|sed "s/,/ /g")

count=0
for np in $NARROWPEAKS;do
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
    awk -F"\t" -v OFS="\t" '$2>0' $np > $TMPDIR/$(basename $np) # avoiding errors like:
# Error: Invalid record in file NHP10_TP4_CD14posDRpos.macs2.narrowPeak.tmp1.bed. Record is
# chr10   -32     468     NHP10_TP4_CD14posDRpos.macs2_peak_16381 13885   .
# Error in read.table(file = file, header = header, sep = sep, quote = quote,  :
#   no lines available in input
# Calls: read.csv -> read.table
# Execution halted
    Rscript ${SCRIPTSFOLDER}/narrowPeak_to_fixed_width_peakSet.R \
    -i $TMPDIR/$(basename $np) \
    -o $outnp \
    -w $WIDTH \
    -t $TMPDIR
done

Rscript ${SCRIPTSFOLDER}/fixed_width_peakSets_to_consensus_peakSet.R \
    -i $OUTNPS \
    -p $REPLICATENAMES \
    -o $CONSENSUSNP \
    -t $TMPDIR \
    -m $CONSENSUSMINREP \
    -c $CONSENSUSMINSPM

Rscript ${SCRIPTSFOLDER}/narrowPeak_normalize_pvalues.R \
    -i $CONSENSUSNP \
    -o $CONSENSUSRENORMNP