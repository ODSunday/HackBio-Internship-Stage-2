#!/bin/bash

# Set the variables
TRIMMED_OUTPUT_DIR="trimmed_reads"         # Directory containing the trimmed reads.
TRIMMED_QC_DIR="qc_trim"                   # Directory for trimmed reads QC reports.
TRIMMED_MULTIQC_DIR="multiqc_trim"         # Directory for the aggregated QC reports.

# Create output directories
mkdir -p "$TRIMMED_QC_DIR"
mkdir -p "$TRIMMED_MULTIQC_DIR"


echo "Performing quality check on the trimmed data..."

fastqc "$TRIMMED_OUTPUT_DIR"/*.trimmed.fastq.gz -o "$TRIMMED_QC_DIR" --quiet
    
echo "Quality check successfully done!"

# Aggregate QC reports with `multiqc`
echo "Aggregating QC reports..."

multiqc "$TRIMMED_QC_DIR" -o "$TRIMMED_MULTIQC_DIR" --quiet

echo "Reports saved to $TRIMMED_MULTIQC_DIR"