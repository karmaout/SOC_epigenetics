#!/usr/bin/env bash
set -euo pipefail

JOB_NAME="mergePseudo"
QUEUE="normal"
NCORES=4
MEM_MB=64000
WALL="12:00"

REPO="/gpfs/home/tqin"
LOG_DIR="${REPO}/logs"
SCRIPT="${REPO}/02_merge_pseudobulk_by_group.sh"

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

       echo WORKDIR=\$(pwd);
       echo SUBSET_COUNT=\$(ls subset_bams/*.bam 2>/dev/null | wc -l);

       bash '$SCRIPT'
     "
