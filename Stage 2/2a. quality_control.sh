#!/bin/bash

# Set the variables
DATA_DIR="raw_data"                        # Directory containing the raw data.
OUTPUT_DIR="qc_reports"                    # Directory for the QC reports.
MULTIQC_DIR="multiqc_reports"              # Directory for the aggregated QC reports.

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$MULTIQC_DIR"

# Quality Control with `fastqc`

echo "Running fastqc..."

for SAMPLE in "$DATA_DIR"/*.fastq.gz; do
    fastqc "$SAMPLE" -o "$OUTPUT_DIR"
done

echo "QC successfully done!"

# Aggregate QC reports with `multiqc`
echo "Aggregating QC reports..."

multiqc "$OUTPUT_DIR" -o "$MULTIQC_DIR"

echo "Reports saved to $MULTIQC_DIR"