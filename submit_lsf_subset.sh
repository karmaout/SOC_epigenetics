#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs

bsub -J "subsetBAM" \
     -n 4 \
     -R "span[hosts=1] rusage[mem=16000]" \
     -M 16000 \
     -W 12:00 \
     -oo logs/subsetBAM.%J.out \
     -eo logs/subsetBAM.%J.err \
     bash -lc 'export THREADS=4; bash scripts/01_subset_bams_by_barcodes.sh'
