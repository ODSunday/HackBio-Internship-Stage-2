#!/bin/bash

# Set variables
REFERENCE_GENOME="./genome/genomeIndex"                      # Reference genome directory.
OUTPUT_DIR="./mapped_reads"                                  # Directory for the STAR output (mapped reads).
READS_DIR="./trimmed_reads"                                  # Directory containing trimmed reads.

# Create output directory
mkdir -p $OUTPUT_DIR

# Create sample array
SAMPLES=(
    "SRR20959676"
    "SRR20959677"
    "SRR20959678"
    "SRR20959679"
    "SRR20959680"
    "SRR20959681"
    "SRR20959682"
)

# Loop through each sample and run STAR
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: $sample"

    # Set input files
    R1="${READS_DIR}/${sample}_1.trimmed.fastq" 
    R2="${READS_DIR}/${sample}_2.trimmed.fastq"

    # Set output prefix
    OUTPUT_PREFIX="${OUTPUT_DIR}/${sample}_"

    echo "Reference Genome: $REFERENCE_GENOME"
    echo "Read 1: $R1"
    echo "Read 2: $R2"
    echo "Output Prefix: $OUTPUT_PREFIX"

    # Run STAR
    STAR --genomeDir $REFERENCE_GENOME --readFilesIn $R1 $R2 --outFileNamePrefix $OUTPUT_PREFIX --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 --outFilterMismatchNmax 2
         
    echo "Finished processing sample: $sample"
done

echo "All samples processed."

# --outSAMtype BAM SortedByCoordinate: Sorts output BAM files by coordinate.
# --outSAMattributes All: Includes all SAM attributes in the output BAM files.
# --outFilterMultimapNmax 1: Keeps uniquely mapped reads only.
# --outFilterMismatchNmax 2: Allows a maximum of 2 mismatches.
