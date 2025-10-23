#!/bin/bash

# Set the variables
DATA_DIR="raw_data"                        # Directory containing the raw data.
TRIMMED_OUTPUT_DIR="trimmed_reads"         # Directory for the trimmed reads.

# Create the output directory
mkdir -p "$TRIMMED_OUTPUT_DIR"

echo "Trimming reads with fastp..."

for FILE in "$DATA_DIR"/*.fastq.gz; do
    base=$(basename "$FILE" .fastq.gz)
    fastp -i "$DATA_DIR/${base}.fastq.gz" \
    -o "$TRIMMED_OUTPUT_DIR/${base}_trimmed.fastq" \
    --html "$TRIMMED_OUTPUT_DIR/${base}_fastp_report.html" --json "$TRIMMED_OUTPUT_DIR/${base}_fastp_report.json"
    --overrepresented_threshold 0.01 \
    --cut_mean_quality 20 --cut_front --cut_tail \
    --qualified_quality_phred 20 --length_required 30 \
    --detect_adapter_for_pe 
done

echo "Results saved into $TRIMMED_OUTPUT_DIR." 