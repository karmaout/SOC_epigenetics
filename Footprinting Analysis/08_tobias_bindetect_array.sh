#!/usr/bin/env bash
set -euo pipefail

TABLE="/gpfs/home/tqin/tobias/bindetect_comparisons.tsv"
FOOT_BASE="/gpfs/home/tqin/tobias/footprints"
OUT_BASE="/gpfs/home/tqin/tobias/bindetect"

MOTIFS="/gpfs/home/tqin/tobias/motifs/JASPAR_vertebrates.meme"
GENOME="/gpfs/home/tqin/mRatBN7/fasta/genome.fa"

THREADS="${THREADS:-4}"
IDX="${LSB_JOBINDEX:?LSB_JOBINDEX not set}"

# Get the IDX-th data row (skip header)
line="$(tail -n +2 "$TABLE" | sed -n "${IDX}p")"
if [ -z "$line" ]; then
  echo "ERROR: No row for index $IDX in $TABLE"
  exit 1
fi

# Parse TSV safely
IFS=$'\t' read -r name prefix cond1 cond2 peaks <<< "$line"

if [ -z "${name:-}" ] || [ -z "${prefix:-}" ] || [ -z "${cond1:-}" ] || [ -z "${cond2:-}" ] || [ -z "${peaks:-}" ]; then
  echo "ERROR: Failed to parse TSV row: $line"
  exit 1
fi

sample1="${prefix}${cond1}"
sample2="${prefix}${cond2}"

sig1="${FOOT_BASE}/${sample1}/footprints.bw"
sig2="${FOOT_BASE}/${sample2}/footprints.bw"
outdir="${OUT_BASE}/${name}"

# Checks
for f in "$sig1" "$sig2" "$peaks" "$MOTIFS" "$GENOME" "${GENOME}.fai"; do
  if [ ! -f "$f" ]; then
    echo "ERROR: Missing required file: $f"
    exit 1
  fi
done

mkdir -p "$outdir"

echo "==== TOBIAS BINDetect ===="
echo "Host   : $(hostname)"
echo "Index  : $IDX"
echo "Name   : $name"
echo "Peaks  : $peaks"
echo "Motifs : $MOTIFS"
echo "Genome : $GENOME"
echo "Cond1  : $cond1 -> $sig1"
echo "Cond2  : $cond2 -> $sig2"
echo "Outdir : $outdir"
echo "Cores  : $THREADS"
echo "=========================="

TOBIAS BINDetect \
  --motifs "$MOTIFS" \
  --signals "$sig1" "$sig2" \
  --genome "$GENOME" \
  --peaks "$peaks" \
  --cond_names "$cond1" "$cond2" \
  --outdir "$outdir" \
  --cores "$THREADS"

echo "DONE: $name"
