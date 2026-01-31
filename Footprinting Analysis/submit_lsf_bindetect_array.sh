#!/usr/bin/env bash
set -euo pipefail

REPO="/gpfs/home/tqin"
TABLE="${REPO}/tobias/bindetect_comparisons.tsv"
LOG_DIR="${REPO}/tobias/logs"
SCRIPT="${REPO}/08_tobias_bindetect_array.sh"

mkdir -p "$LOG_DIR"

N=$(($(wc -l < "$TABLE") - 1))
if [ "$N" -le 0 ]; then
  echo "ERROR: No comparisons in $TABLE"
  exit 1
fi

QUEUE="normal"
NCORES=4
MEM_MB=64000
WALL="12:00"

bsub -q "$QUEUE" \
     -J "BINDetect[1-${N}]" \
     -n "$NCORES" \
     -R "span[hosts=1] rusage[mem=${MEM_MB}]" \
     -M "$MEM_MB" \
     -W "$WALL" \
     -oo "${LOG_DIR}/BINDetect.%J.%I.out" \
     -eo "${LOG_DIR}/BINDetect.%J.%I.err" \
     bash -lc "cd '$REPO';
       source ~/miniconda3/etc/profile.d/conda.sh;
       conda activate tobias_ppc;
       export THREADS=${NCORES};
       bash '$SCRIPT'
     "
