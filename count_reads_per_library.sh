#!/bin/bash

# Output file
OUTPUT="library_read_counts.tsv"
echo -e "library_id\tread_count" > $OUTPUT

# Define your libraries
LIBRARIES=("02BLA" "29BLA" "79BLA" "81BLA" "82BLA" "85BLA" "89BLA" "90BLA" "91BLA" "93BLA" "94BLA" "95BLA"
           "02PC" "29PC" "79PC" "81PC" "82PC" "85PC" "89PC" "90PC" "91PC" "93PC" "94PC" "95PC")

# Base directory where your BAM files are located
BAM_DIR="/mnt/DATA/Tian/SOC"

for LIB in "${LIBRARIES[@]}"; do
  BAM_FILE="$BAM_DIR/$LIB/outs/atac_possorted_bam.bam"
  if [[ -f "$BAM_FILE" ]]; then
    COUNT=$(samtools idxstats "$BAM_FILE" | awk '{s+=$3} END {print s}')
    echo -e "${LIB}\t${COUNT}" >> $OUTPUT
    echo "Processed $LIB: $COUNT reads"
  else
    echo "Warning: BAM not found for $LIB at $BAM_FILE"
  fi
done

echo "✅ All counts saved to $OUTPUT"
