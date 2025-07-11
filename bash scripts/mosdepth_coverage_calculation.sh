# Path setup
BAMDIR="merged_bams"
BED="bins_50kb_3col.bed"

# Output prefix list
GROUPS=("P_E" "P_R" "U_E" "U_R")

for GROUP in "${GROUPS[@]}"; do
  echo "📊 Running mosdepth for $GROUP..."
  mosdepth --by "$BED" "$GROUP" "${BAMDIR}/${GROUP}_merged.sorted.bam"
  zcat "${GROUP}.regions.bed.gz" > "${GROUP}_coverage_50kb.tsv"
done
