#!/usr/bin/env bash
set -euo pipefail

# =========================
# LSF job configuration
# =========================
JOB_NAME="subsetBAM"
QUEUE="normal"
NCORES=4
MEM_MB=64000      # 64 GB (safe for large ATAC BAMs)
WALL="24:00"

# =========================
# Paths
# =========================
LOG_DIR="logs"
PIPELINE_SCRIPT="01_subset_bams_by_barcodes.sh"

mkdir -p "$LOG_DIR"

echo "Submitting LSF job:"
echo "  Job name : $JOB_NAME"
echo "  Cores    : $NCORES"
echo "  Memory   : ${MEM_MB} MB"
echo "  Walltime : $WALL"
echo "  Script   : $PIPELINE_SCRIPT"
echo

bsub -q "$QUEUE" \
     -J "$JOB_NAME" \
     -n "$NCORES" \
     -R "span[hosts=1] rusage[mem=${MEM_MB}]" \
     -M "$MEM_MB" \
     -W "$WALL" \
     -oo "${LOG_DIR}/${JOB_NAME}.%J.out" \
     -eo "${LOG_DIR}/${JOB_NAME}.%J.err" \
     bash -lc "
       source ~/miniconda3/etc/profile.d/conda.sh
       conda activate tobias_ppc
       export THREADS=${NCORES}
       bash ${PIPELINE_SCRIPT}
     "
