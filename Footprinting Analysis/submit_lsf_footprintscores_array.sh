#!/usr/bin/env bash
set -euo pipefail

REPO="/gpfs/home/tqin"
BAM_LIST="${REPO}/tobias/bam_list_atacorrect.txt"
LOG_DIR="${REPO}/tobias/logs"
SCRIPT="${REPO}/07_tobias_footprintscores_array.sh"

mkdir -p "$LOG_DIR"

if [ ! -f "$BAM_LIST" ]; then
  echo "ERROR: Missing $BAM_LIST"
  echo "Create it with: ls -1 ${REPO}/pseudobulk_bams/*.bam | sort > $BAM_LIST"
  exit 1
fi

N=$(wc -l < "$BAM_LIST")
if [ "$N" -eq 0 ]; then
  echo "ERROR: $BAM_LIST is empty"
  exit 1
fi

QUEUE="normal"
NCORES=4
MEM_MB=64000
WALL="12:00"

bsub -q "$QUEUE" \
     -J "Footprints[1-${N}]" \
     -n "$NCORES" \
     -R "span[hosts=1] rusage[mem=${MEM_MB}]" \
     -M "$MEM_MB" \
     -W "$WALL" \
     -oo "${LOG_DIR}/Footprints.%J.%I.out" \
     -eo "${LOG_DIR}/Footprints.%J.%I.err" \
     bash -lc "cd '$REPO';
       source ~/miniconda3/etc/profile.d/conda.sh;
       conda activate tobias_ppc;
       export THREADS=${NCORES};
       bash '$SCRIPT'
     "
