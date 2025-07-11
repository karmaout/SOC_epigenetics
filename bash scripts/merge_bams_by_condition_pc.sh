#!/bin/sh

# Output directory for PC (absolute path)
OUTDIR="$(pwd)/merged_bams_PC"
mkdir -p "$OUTDIR"

# Create BAM list files by behavior group
echo "/mnt/DATA/Tian/SOC/02PC/outs/atac_possorted_bam.bam" >  "$OUTDIR/P_R.txt"
echo "/mnt/DATA/Tian/SOC/82PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_R.txt"
echo "/mnt/DATA/Tian/SOC/91PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_R.txt"

echo "/mnt/DATA/Tian/SOC/29PC/outs/atac_possorted_bam.bam" >  "$OUTDIR/P_E.txt"
echo "/mnt/DATA/Tian/SOC/79PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_E.txt"
echo "/mnt/DATA/Tian/SOC/81PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_E.txt"

echo "/mnt/DATA/Tian/SOC/85PC/outs/atac_possorted_bam.bam" >  "$OUTDIR/U_E.txt"
echo "/mnt/DATA/Tian/SOC/89PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_E.txt"
echo "/mnt/DATA/Tian/SOC/94PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_E.txt"

echo "/mnt/DATA/Tian/SOC/90PC/outs/atac_possorted_bam.bam" >  "$OUTDIR/U_R.txt"
echo "/mnt/DATA/Tian/SOC/93PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_R.txt"
echo "/mnt/DATA/Tian/SOC/95PC/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_R.txt"

# Merge, sort, and index each group
for GROUP in P_E P_R U_E U_R; do
  echo "🔄 Merging BAMs for group $GROUP..."
  samtools merge -@ 4 -b "$OUTDIR/${GROUP}.txt" "$OUTDIR/${GROUP}_merged.bam"
  samtools sort -@ 4 -o "$OUTDIR/${GROUP}_merged.sorted.bam" "$OUTDIR/${GROUP}_merged.bam"
  samtools index "$OUTDIR/${GROUP}_merged.sorted.bam"
  rm "$OUTDIR/${GROUP}_merged.bam"
done

echo "✅ All BAMs for PC merged and indexed."
echo "📂 Output directory: $OUTDIR"

