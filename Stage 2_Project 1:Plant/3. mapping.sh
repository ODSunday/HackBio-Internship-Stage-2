#!/bin/bash

# Set variables
REFERENCE_GENOME="./genome/genomeIndex"                      # Reference genome directory.
OUTPUT_DIR="./mapped_reads"                                  # Directory for the STAR output (mapped reads).
READS_DIR="./trimmed_reads"                                  # Directory containing trimmed reads.

# Create output directory
mkdir -p $OUTPUT_DIR

# Create sample array
SAMPLES=(
    "SRR12808527_c1"
    "SRR12808528_c2"
    "SRR12808529_c3"
    "SRR12808497_t1"
    "SRR12808498_t2"
    "SRR12808499_t3"
)

# Loop through each sample and run STAR
for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: $sample"

    # Set input files
    infile="${READS_DIR}/${sample}_trimmed.fastq" 
    
    # Set output prefix
    OUTPUT_PREFIX="${OUTPUT_DIR}/${sample}_"

    echo "Reference Genome: $REFERENCE_GENOME"
    echo "Input file: $infile"
    echo "Output Prefix: $OUTPUT_PREFIX"

    # Run STAR
    STAR --genomeDir $REFERENCE_GENOME --readFilesIn $infile --outFileNamePrefix $OUTPUT_PREFIX --outSAMtype BAM SortedByCoordinate --outSAMattributes All --outFilterMultimapNmax 1 --outFilterMismatchNmax 2
         
    echo "Finished processing sample: $sample"
done

echo "All samples processed."

# --outSAMtype BAM SortedByCoordinate: Sorts output BAM files by coordinate.
# --outSAMattributes All: Includes all SAM attributes in the output BAM files.
# --outFilterMultimapNmax 1: Keeps uniquely mapped reads only.
# --outFilterMismatchNmax 2: Allows a maximum of 2 mismatches.