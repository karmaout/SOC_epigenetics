#!/usr/bin/env bash
set -euo pipefail

REPO="/gpfs/home/tqin"
BAM_DIR="${REPO}/pseudobulk_bams"
LIST_DIR="${REPO}/tobias"
BAM_LIST="${LIST_DIR}/bam_list_atacorrect.txt"
LOG_DIR="${LIST_DIR}/logs"

SCRIPT="${REPO}/06_tobias_atacorrect_array.sh"

mkdir -p "$LIST_DIR" "$LOG_DIR"

# Build BAM list (sorted for reproducibility)
ls -1 ${BAM_DIR}/*.bam | sort > "$BAM_LIST"
N=$(wc -l < "$BAM_LIST")

if [ "$N" -eq 0 ]; then
  echo "ERROR: No BAMs found in ${BAM_DIR}"
  exit 1
fi

echo "Found ${N} BAMs for TOBIAS ATACorrect."
echo "BAM list: $BAM_LIST"
head -n 20 "$BAM_LIST"

# ---- LSF resources ----
QUEUE="normal"
NCORES=4
MEM_MB=64000
WALL="12:00"

bsub -q "$QUEUE" \
     -J "ATACorrect[1-${N}]" \
     -n "$NCORES" \
     -R "span[hosts=1] rusage[mem=${MEM_MB}]" \
     -M "$MEM_MB" \
     -W "$WALL" \
     -oo "${LOG_DIR}/ATACorrect.%J.%I.out" \
     -eo "${LOG_DIR}/ATACorrect.%J.%I.err" \
     bash -lc "cd '$REPO';
       source ~/miniconda3/etc/profile.d/conda.sh;
       conda activate tobias_ppc;
       export THREADS=${NCORES};
       bash '$SCRIPT'
     "
