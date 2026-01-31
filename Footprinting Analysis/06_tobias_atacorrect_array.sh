#!/usr/bin/env bash
set -euo pipefail

# ---------- Paths ----------
BAM_LIST="/gpfs/home/tqin/tobias/bam_list_atacorrect.txt"
GENOME="/gpfs/home/tqin/mRatBN7/fasta/genome.fa"
OUT_BASE="/gpfs/home/tqin/tobias/atacorrect"

BLA_PEAKS="/gpfs/home/tqin/peaks/BLA_VGLUT1_union_peaks.bed"
PC_PEAKS="/gpfs/home/tqin/peaks/PC_VGLUT1_union_peaks.bed"

THREADS="${THREADS:-4}"

# ---------- LSF array index ----------
IDX="${LSB_JOBINDEX:?LSB_JOBINDEX not set (run as LSF array job)}"
bam="$(sed -n "${IDX}p" "$BAM_LIST")"

if [ -z "$bam" ] || [ ! -f "$bam" ]; then
  echo "ERROR: BAM not found for index ${IDX}: '${bam}'"
  exit 1
fi

sample="$(basename "$bam" .bam)"

# Choose peaks by prefix
peaks=""
if [[ "$sample" == BLA_VGLUT1__* ]]; then
  peaks="$BLA_PEAKS"
elif [[ "$sample" == PC_VGLUT1__* ]]; then
  peaks="$PC_PEAKS"
else
  echo "ERROR: Unknown sample prefix for '${sample}' (expected BLA_VGLUT1__* or PC_VGLUT1__*)"
  exit 1
fi

# Sanity checks
for f in "$GENOME" "$GENOME.fai" "$peaks"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Required file missing: $f"
    exit 1
  fi
done

outdir="${OUT_BASE}/${sample}"
mkdir -p "$outdir"

echo "==== TOBIAS ATACorrect ===="
echo "Host   : $(hostname)"
echo "Index  : ${IDX}"
echo "Sample : ${sample}"
echo "BAM    : ${bam}"
echo "Peaks  : ${peaks}"
echo "Genome : ${GENOME}"
echo "Outdir : ${outdir}"
echo "Cores  : ${THREADS}"
echo "==========================="

TOBIAS ATACorrect \
  --bam "$bam" \
  --genome "$GENOME" \
  --peaks "$peaks" \
  --outdir "$outdir" \
  --cores "$THREADS"

echo "DONE: ${sample}"
