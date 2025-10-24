#!/bin/bash

# Set the path to the GTF file
ANNOTATION_FILE="./genome/a_thaliana.gff3"                          # Path to the annotation file.
OUTPUT_DIR="./counts"                                               # Directory for featureCounts output.
MAPPED_READS_DIR="./mapped_reads"                                   # Directory containing the mapped reads.

# Create output directory
mkdir -p $OUTPUT_DIR

# Collect all BAM files from the `mapped_reads` directory
BAM_FILES=(${MAPPED_READS_DIR}/*.bam)

# Check if BAM files are present
if [ ${#BAM_FILES[@]} -eq 0 ]; then
    echo "No BAM files found in $MAPPED_READS_DIR"
    exit 1
else 
    echo "BAM files present. Proceed with counting..."
fi

# Output file for counts
OUTPUT_COUNT_FILE="$OUTPUT_DIR/counts.txt"

# Run featureCounts
echo "Running featureCounts..."
featureCounts -a $ANNOTATION_FILE -o $OUTPUT_COUNT_FILE -t gene -g gene_id -C -O "${BAM_FILES[@]}"

# Check if featureCounts ran successfully
if [ $? -eq 0 ]; then
    echo "Counting completed successfully. Output saved to $OUTPUT_COUNT_FILE"
else
    echo "Error occurred during featureCounts execution."
    exit 1
fi 

# -a: annotation file (GFF3)
# -o: output file
# -t: feature type to count (gene)
# -g: attribute to group features (gene_id)
# -C: counts only reads that are uniquely mapped
# -O: counts reads overlapping multiple features only once - avoiding double-counting
