#!/usr/bin/env bash
set -euo pipefail

# ---- Inputs ----
BAM_LIST="/gpfs/home/tqin/peaks/bam_list.txt"
OUTDIR="/gpfs/home/tqin/peaks/macs2"
QVAL="0.01"

# For ATAC, these are standard choices:
SHIFT="-100"
EXTSIZE="200"

# MACS2 genome size:
# -g mm is commonly used (~1.87e9). For rat, MACS2 does not have a perfect builtin.
# If you prefer explicit, use: -g 2.5e9 (rat rough) or your exact mappable genome size.
GENOME_SIZE="2.5e9"

# ---- LSF array index ----
IDX="${LSB_JOBINDEX:?LSB_JOBINDEX not set (run as LSF array job)}"

bam="$(sed -n "${IDX}p" "$BAM_LIST")"
if [ -z "$bam" ]; then
  echo "ERROR: No BAM found for index $IDX in $BAM_LIST"
  exit 1
fi

sample="$(basename "$bam" .bam)"

mkdir -p "$OUTDIR"

echo "==== MACS2 peak calling ===="
echo "Host   : $(hostname)"
echo "Index  : $IDX"
echo "BAM    : $bam"
echo "Sample : $sample"
echo "Outdir : $OUTDIR"
echo "============================"

# Call peaks
macs2 callpeak \
  -t "$bam" \
  -f BAM \
  -g "$GENOME_SIZE" \
  --nomodel \
  --shift "$SHIFT" \
  --extsize "$EXTSIZE" \
  -q "$QVAL" \
  -n "$sample" \
  --outdir "$OUTDIR"

echo "DONE: $sample"
