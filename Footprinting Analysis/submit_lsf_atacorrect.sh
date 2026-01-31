#!/usr/bin/env bash
set -euo pipefail

JOB_NAME="ATACorrect"
QUEUE="normal"
NCORES=4
MEM_MB=64000
WALL="08:00"

REPO="/gpfs/home/tqin"
LOG_DIR="${REPO}/tobias/logs"
SCRIPT="${REPO}/03_tobias_atacorrect.sh"

mkdir -p "$LOG_DIR"

bsub -q "$QUEUE" \
     -J "$JOB_NAME" \
     -n "$NCORES" \
     -R "span[hosts=1] rusage[mem=${MEM_MB}]" \
     -M "$MEM_MB" \
     -W "$WALL" \
     -oo "${LOG_DIR}/${JOB_NAME}.%J.out" \
     -eo "${LOG_DIR}/${JOB_NAME}.%J.err" \
     bash -lc "set -e;
       cd '$REPO';
       source ~/miniconda3/etc/profile.d/conda.sh;
       conda activate tobias_ppc;
       export THREADS=$NCORES;
       bash '$SCRIPT'
     "
