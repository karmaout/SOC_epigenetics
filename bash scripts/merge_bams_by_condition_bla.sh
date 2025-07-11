
#!/bin/sh

# Output directory (absolute path)
OUTDIR="$(pwd)/merged_bams"
mkdir -p "$OUTDIR"

# Create BAM list files by behavior group
echo "/mnt/DATA/Tian/SOC/02BLA/outs/atac_possorted_bam.bam" >  "$OUTDIR/P_R.txt"
echo "/mnt/DATA/Tian/SOC/82BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_R.txt"
echo "/mnt/DATA/Tian/SOC/91BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_R.txt"

echo "/mnt/DATA/Tian/SOC/29BLA/outs/atac_possorted_bam.bam" >  "$OUTDIR/P_E.txt"
echo "/mnt/DATA/Tian/SOC/79BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_E.txt"
echo "/mnt/DATA/Tian/SOC/81BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/P_E.txt"

echo "/mnt/DATA/Tian/SOC/85BLA/outs/atac_possorted_bam.bam" >  "$OUTDIR/U_E.txt"
echo "/mnt/DATA/Tian/SOC/89BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_E.txt"
echo "/mnt/DATA/Tian/SOC/94BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_E.txt"

echo "/mnt/DATA/Tian/SOC/90BLA/outs/atac_possorted_bam.bam" >  "$OUTDIR/U_R.txt"
echo "/mnt/DATA/Tian/SOC/93BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_R.txt"
echo "/mnt/DATA/Tian/SOC/95BLA/outs/atac_possorted_bam.bam" >> "$OUTDIR/U_R.txt"

# Merge, sort, and index each group
for GROUP in P_E P_R U_E U_R; do
  echo "🔄 Merging BAMs for group $GROUP..."
  samtools merge -@ 4 -b "$OUTDIR/${GROUP}.txt" "$OUTDIR/${GROUP}_merged.bam"
  samtools sort -@ 4 -o "$OUTDIR/${GROUP}_merged.sorted.bam" "$OUTDIR/${GROUP}_merged.bam"
  samtools index "$OUTDIR/${GROUP}_merged.sorted.bam"
  rm "$OUTDIR/${GROUP}_merged.bam"
done

echo "✅ All BAMs merged and indexed."
echo "📂 Output directory: $OUTDIR"
