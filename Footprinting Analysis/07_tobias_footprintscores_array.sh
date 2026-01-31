#!/usr/bin/env bash
set -euo pipefail

# ----- Inputs -----
BAM_LIST="/gpfs/home/tqin/tobias/bam_list_atacorrect.txt"
CORR_BASE="/gpfs/home/tqin/tobias/atacorrect"
OUT_BASE="/gpfs/home/tqin/tobias/footprints"

# Union peaks (region-specific)
BLA_PEAKS="/gpfs/home/tqin/peaks/BLA_VGLUT1_union_peaks.bed"
PC_PEAKS="/gpfs/home/tqin/peaks/PC_VGLUT1_union_peaks.bed"

THREADS="${THREADS:-4}"
IDX="${LSB_JOBINDEX:?LSB_JOBINDEX not set (run as LSF array job)}"

# ----- Select BAM/sample by array index -----
bam="$(sed -n "${IDX}p" "$BAM_LIST")"
if [ -z "$bam" ] || [ ! -f "$bam" ]; then
  echo "ERROR: BAM not found for index ${IDX}: '${bam}'"
  exit 1
fi

sample="$(basename "$bam" .bam)"

# ----- Choose peaks by sample prefix -----
if [[ "$sample" == BLA_VGLUT1__* ]]; then
  peaks="$BLA_PEAKS"
elif [[ "$sample" == PC_VGLUT1__* ]]; then
  peaks="$PC_PEAKS"
else
  echo "ERROR: Unknown sample prefix: $sample"
  exit 1
fi

# ----- Paths -----
corr_bw="${CORR_BASE}/${sample}/${sample}_corrected.bw"
outdir="${OUT_BASE}/${sample}"
mkdir -p "$outdir"

# ----- Checks -----
for f in "$corr_bw" "$peaks"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

echo "==== TOBIAS FootprintScores ===="
echo "Host   : $(hostname)"
echo "Index  : $IDX"
echo "Sample : $sample"
echo "BW     : $corr_bw"
echo "Peaks  : $peaks"
echo "Outdir : $outdir"
echo "Cores  : $THREADS"
echo "==============================="

TOBIAS FootprintScores \
  --signal "$corr_bw" \
  --regions "$peaks" \
  --output "${outdir}/footprints.bw" \
  --cores "$THREADS"

echo "DONE: $sample"
