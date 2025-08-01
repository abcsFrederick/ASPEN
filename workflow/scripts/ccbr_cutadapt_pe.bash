#!/bin/bash

# . /opt2/conda/etc/profile.d/conda.sh
# conda activate python3

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="Trim PE reads using cutadapt"      # this is optional
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--infastq1',required=True, help='input R1 fastq.gz file')
parser.add_argument('--infastq2',required=True, help='input R2 fastq.gz file')
parser.add_argument('--samplename',required=True, help='samplename')
parser.add_argument('--threads',required=True, help='number of threads')
parser.add_argument('--outfastq1',required=True, help='output R1 fastq.gz file')
parser.add_argument('--outfastq2',required=True, help='output R2 fastq.gz file')
# parser.add_argument('--sampleinfo',required=True, help='nreads file')
EOF

infastq1=$INFASTQ1
infastq2=$INFASTQ2
samplename=$SAMPLENAME
outfastq1=$OUTFASTQ1
outfastq2=$OUTFASTQ2
ncpus=$THREADS

cutadapt \
--pair-filter=any \
--nextseq-trim=2 \
--trim-n \
-n 5 -O 5 \
-q 10,10 \
-m 35:35 \
-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
-B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
-j $ncpus \
-o $outfastq1 \
-p $outfastq2 \
$infastq1 $infastq2

# nreads=`zcat $infastq1|wc -l`
# nreadstrim=`zcat $outfastq1|wc -l`
# echo "$nreads $nreadstrim"|awk '{printf("%d\tInput Nreads\n%d\tAfter trimming\n",$1/2,$2/2)}' > ${SAMPLEINFO}

# conda deactivate
