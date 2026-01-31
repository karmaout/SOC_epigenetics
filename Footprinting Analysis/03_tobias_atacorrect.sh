#!/usr/bin/env bash
set -euo pipefail

GENOME="/gpfs/home/tqin/mRatBN7/fasta/genome.fa"
IN_DIR="/gpfs/home/tqin/pseudobulk_bams"
OUT_BASE="/gpfs/home/tqin/tobias/atacorrect"
THREADS="${THREADS:-4}"

mkdir -p "$OUT_BASE"

for bam in ${IN_DIR}/*.bam; do
  sample=$(basename "$bam" .bam)
  outdir="${OUT_BASE}/${sample}"

  echo
  echo "===== ATACorrect: ${sample} ====="
  echo "BAM    : $bam"
  echo "GENOME : $GENOME"
  echo "OUTDIR : $outdir"
  echo "CORES  : $THREADS"

  TOBIAS ATACorrect \
    --bam "$bam" \
    --genome "$GENOME" \
    --outdir "$outdir" \
    --cores "$THREADS"
done
