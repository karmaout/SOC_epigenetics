#!/usr/bin/env bash
set -euo pipefail

PEAK_DIR="/gpfs/home/tqin/peaks/macs2"
OUT_DIR="/gpfs/home/tqin/peaks"

mkdir -p "$OUT_DIR"

merge_bed() {
  # stdin: sorted BED3 (chr, start, end)
  # stdout: merged BED3
  awk '
    BEGIN{OFS="\t"}
    NR==1 {c=$1; s=$2; e=$3; next}
    $1==c && $2<=e { if($3>e) e=$3; next }
    {print c,s,e; c=$1; s=$2; e=$3}
    END{print c,s,e}
  '
}

make_union_for_prefix() {
  local prefix="$1"
  local outbed="$2"

  files=( "${PEAK_DIR}/${prefix}"*_peaks.narrowPeak )

  if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: No narrowPeak files found for prefix '${prefix}' in ${PEAK_DIR}"
    exit 1
  fi

  echo "Building union peaks for: ${prefix}"
  echo "  Inputs: ${#files[@]} files"
  echo "  Output: ${outbed}"

  # narrowPeak columns: chrom start end ...
  # Convert to BED3, sort, merge overlaps/touching.
  cat "${files[@]}" \
    | awk 'BEGIN{OFS="\t"}{print $1,$2,$3}' \
    | sort -k1,1 -k2,2n \
    | merge_bed \
    > "${outbed}"

  # Basic stats
  n=$(wc -l < "${outbed}" | awk '{print $1}')
  echo "  Union peaks: ${n}"
  echo
}

# ----- Make union peak sets -----
make_union_for_prefix "BLA_VGLUT1__" "${OUT_DIR}/BLA_VGLUT1_union_peaks.bed"
make_union_for_prefix "PC_VGLUT1__"  "${OUT_DIR}/PC_VGLUT1_union_peaks.bed"

# Optional: global union across both (useful for a single mask)
cat "${OUT_DIR}/BLA_VGLUT1_union_peaks.bed" "${OUT_DIR}/PC_VGLUT1_union_peaks.bed" \
  | sort -k1,1 -k2,2n \
  | merge_bed \
  > "${OUT_DIR}/VGLUT1_BLA_PC_union_peaks.bed"

echo "Global union written to: ${OUT_DIR}/VGLUT1_BLA_PC_union_peaks.bed"
echo "Done: $(date)"
