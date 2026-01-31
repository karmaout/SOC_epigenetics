#!/usr/bin/env bash
set -euo pipefail

REPO="/gpfs/home/tqin"
BAM_DIR="${REPO}/pseudobulk_bams"
LIST_DIR="${REPO}/peaks"
BAM_LIST="${LIST_DIR}/bam_list.txt"
LOG_DIR="${LIST_DIR}/logs"
SCRIPT="${REPO}/04_callpeaks_macs2.sh"

mkdir -p "$LIST_DIR" "$LOG_DIR"

# Create BAM list
ls -1 ${BAM_DIR}/*.bam > "$BAM_LIST"
N=$(wc -l < "$BAM_LIST")

if [ "$N" -eq 0 ]; then
  echo "ERROR: No BAMs found in ${BAM_DIR}"
  exit 1
fi

echo "Found ${N} BAMs for peak calling."
echo "BAM list: $BAM_LIST"
head -n 5 "$BAM_LIST"

# ---- LSF resources ----
QUEUE="normal"
NCORES=4
MEM_MB=64000
WALL="12:00"

# Submit array job
bsub -q "$QUEUE" \
     -J "callPeaks[1-${N}]" \
     -n "$NCORES" \
     -R "span[hosts=1] rusage[mem=${MEM_MB}]" \
     -M "$MEM_MB" \
     -W "$WALL" \
     -oo "${LOG_DIR}/callPeaks.%J.%I.out" \
     -eo "${LOG_DIR}/callPeaks.%J.%I.err" \
     bash -lc "cd '$REPO';
       source ~/miniconda3/etc/profile.d/conda.sh;
       conda activate tobias_ppc;
       which macs2;
       bash '$SCRIPT'
     "
