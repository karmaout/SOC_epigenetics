#!/usr/bin/env bash
set -euo pipefail

SUBSET_DIR="${SUBSET_DIR:-subset_bams}"
OUTDIR="${OUTDIR:-pseudobulk_bams}"
THREADS="${THREADS:-4}"

PREFIXES=("BLA_VGLUT1" "PC_VGLUT1")
BEHAV_GROUPS=("U_A" "P_A" "U_R" "P_R")


mkdir -p "$OUTDIR"

echo "==== Merge pseudo-bulk BAMs ===="
echo "Start   : $(date)"
echo "Host    : $(hostname)"
echo "THREADS : $THREADS"
echo "SUBSET  : $SUBSET_DIR"
echo "OUTDIR  : $OUTDIR"
echo "================================"

for prefix in "${PREFIXES[@]}"; do
  for grp in "${BEHAV_GROUPS[@]}"; do

    shopt -s nullglob
    files=( "${SUBSET_DIR}/${prefix}__${grp}__"*.bam )
    shopt -u nullglob

    if [ ${#files[@]} -eq 0 ]; then
      echo "SKIP: No files for ${prefix} ${grp}"
      continue
    fi

    out="${OUTDIR}/${prefix}__${grp}.bam"

    echo
    echo "[$(date)] Merging ${prefix} ${grp}"
    echo "  Inputs: ${#files[@]} BAMs"
    echo "  Output: $out"

    samtools merge -f "$out" "${files[@]}"
    samtools sort -@ "$THREADS" -o "$out" "$out"
    samtools index "$out"

    n=$(samtools idxstats "$out" | awk '{s+=$3} END{print s+0}')
    echo "  Merged mapped reads: $n"
  done
done

echo
echo "Done    : $(date)"
