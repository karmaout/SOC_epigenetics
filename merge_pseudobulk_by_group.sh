#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../config/env.sh"

mkdir -p "${PSEUDOBULK_DIR}" "${LOG_DIR}"

echo "Start: $(date)"
echo "Host: $(hostname)"
echo "THREADS=${THREADS}"
echo "SUBSET_DIR=${SUBSET_DIR}"
echo "PSEUDOBULK_DIR=${PSEUDOBULK_DIR}"
echo

for prefix in "${PREFIXES[@]}"; do
  for grp in "${GROUPS[@]}"; do
    pattern="${SUBSET_DIR}/${prefix}__${grp}__"*.bam
    files=( $pattern )

    if [ ${#files[@]} -eq 0 ]; then
      echo "SKIP: No files for ${prefix} ${grp}"
      continue
    fi

    out="${PSEUDOBULK_DIR}/${prefix}__${grp}.bam"
    echo "[$(date)] Merging: ${prefix} ${grp}"
    echo "  -> ${out}"
    echo "  Inputs: ${#files[@]} BAMs"

    samtools merge -f -o "$out" "${files[@]}"
    samtools sort -@ "${THREADS}" -o "$out" "$out"
    samtools index "$out"

    n=$(samtools idxstats "$out" | awk '{s+=$3} END{print s+0}')
    echo "  merged mapped reads: ${n}"
    echo
  done
done

echo "End: $(date)"
