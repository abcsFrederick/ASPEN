#!/bin/bash

set -e -x -o pipefail
ARGPARSE_DESCRIPTION="calculate FLD from dedup bam file" 
source /opt2/argparse.bash || exit 1
argparse "$@" <<EOF || exit 1
parser.add_argument('--dedupbam',required=True, help='dedupbam file')
parser.add_argument('--fldout',required=True, help='output FLD file')
parser.add_argument('--scriptsfolder',required=True, help='folder where the scripts are... generally workflow/scripts')
EOF

python ${SCRIPTSFOLDER}/ccbr_atac_bam2FLD.py -i $DEDUPBAM -o $FLDOUT
