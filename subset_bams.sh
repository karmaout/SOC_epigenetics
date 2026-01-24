#!/usr/bin/env bash
set -euo pipefail

BASE="/gpfs/home/tqin"
BARCODE_DIR="barcode_lists"     # one barcode per line
OUTDIR="subset_bams"

mkdir -p "$OUTDIR"

echo "Start: $(date)"
echo "Running on: $(hostname)"
echo "Conda env: ${CONDA_DEFAULT_ENV:-NA}"

for bcfile in ${BARCODE_DIR}/*.txt; do
  base=$(basename "$bcfile" .txt)   # e.g. BLA_VGLUT1__U_R__95BLA
  lib=${base##*__}
  inbam="${BASE}/${lib}/outs/atac_possorted_bam.bam"
  outbam="${OUTDIR}/${base}.bam"

  if [ ! -f "$inbam" ]; then
    echo "Missing BAM for ${lib}: $inbam"
    continue
  fi

  echo "[$(date)] Subsetting: $base"
  echo "  inbam : $inbam"
  echo "  outbam: $outbam"

  # subset by CB tag using awk hash lookup
  samtools view -h "$inbam" \
  | awk -v BC="$bcfile" '
      BEGIN{
        while((getline line < BC)>0){
          gsub(/\r/,"",line);
          keep["CB:Z:" line]=1
        }
      }
      /^@/ {print; next}
      {
        for(i=12;i<=NF;i++){
          if($i ~ /^CB:Z:/ && keep[$i]) {print; break}
        }
      }' \
  | samtools view -b -o "$outbam" -

  samtools sort -@ "${THREADS:-1}" -o "$outbam" "$outbam"
  samtools index "$outbam"

  # quick read count
  n=$(samtools idxstats "$outbam" | awk '{s+=$3} END{print s+0}')
  echo "  subset mapped reads: $n"
done

echo "End: $(date)"
