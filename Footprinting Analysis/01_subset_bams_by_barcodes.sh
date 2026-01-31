#!/usr/bin/env bash
set -euo pipefail

# =========================
# User-configurable settings
# =========================

# Base directory containing per-library folders like:
# /gpfs/home/tqin/95PC/outs/atac_possorted_bam.bam
BASE="${BASE:-/gpfs/home/tqin}"

# Barcode lists directory (fixed barcodes ending with -1)
BARCODE_DIR="${BARCODE_DIR:-barcode_lists_fixed}"

# Output directory for per-library subset BAMs
OUTDIR="${OUTDIR:-subset_bams}"

# Threads for samtools sort
THREADS="${THREADS:-1}"

# =========================
# Setup
# =========================
mkdir -p "$OUTDIR"

echo "==== Subset BAMs by barcodes ===="
echo "Start   : $(date)"
echo "Host    : $(hostname)"
echo "BASE    : $BASE"
echo "BARCODE : $BARCODE_DIR"
echo "OUTDIR  : $OUTDIR"
echo "THREADS : $THREADS"
echo "================================="

shopt -s nullglob
bc_files=("$BARCODE_DIR"/*.txt)
if [ ${#bc_files[@]} -eq 0 ]; then
  echo "ERROR: No barcode files found in $BARCODE_DIR (expected *.txt)"
  exit 1
fi

# =========================
# Main loop
# =========================
for bcfile in "$BARCODE_DIR"/*.txt; do
  base=$(basename "$bcfile" .txt)          # e.g. BLA_VGLUT1__U_A__85BLA
  lib=${base##*__}                         # e.g. 85BLA
  inbam="${BASE}/${lib}/outs/atac_possorted_bam.bam"
  outbam="${OUTDIR}/${base}.bam"

  if [ ! -f "$inbam" ]; then
    echo "WARN: Missing input BAM for ${lib}: $inbam"
    continue
  fi

  echo
  echo "[$(date)] Processing: $base"
  echo "  Library : $lib"
  echo "  Input   : $inbam"
  echo "  Barcodes: $bcfile"
  echo "  Output  : $outbam"

  # Subset reads where any field equals "CB:Z:<barcode>" for barcodes in bcfile
  # - removes potential Windows CR (\r)
  # - trims trailing whitespace
  samtools view -h "$inbam" \
  | awk -v BC="$bcfile" '
      BEGIN{
        while((getline line < BC)>0){
          gsub(/\r/,"",line);
          gsub(/[[:space:]]+$/,"",line);
          if(line != "") keep["CB:Z:" line]=1;
        }
      }
      /^@/ {print; next}
      {
        for(i=12;i<=NF;i++){
          if($i ~ /^CB:Z:/ && keep[$i]) {print; break}
        }
      }' \
  | samtools view -b -o "$outbam" -

  # Sort + index (required for downstream merging and TOBIAS)
  samtools sort -@ "$THREADS" -o "$outbam" "$outbam"
  samtools index "$outbam"

  # Quick count of mapped reads in output
  n=$(samtools idxstats "$outbam" | awk '{s+=$3} END{print s+0}')
  echo "  Output mapped reads: $n"
done

echo
echo "Done    : $(date)"
