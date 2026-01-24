#!/usr/bin/env bash
set -euo pipefail

mkdir -p logs

bsub -J "mergePseudo" \
     -n 4 \
     -R "span[hosts=1] rusage[mem=16000]" \
     -M 16000 \
     -W 04:00 \
     -oo logs/mergePseudo.%J.out \
     -eo logs/mergePseudo.%J.err \
     bash -lc 'export THREADS=4; bash scripts/02_merge_pseudobulk_by_group.sh'
