#!/bin/bash

INDEX_DIR="./igv"

echo "Indexing in progress"
echo "--------------------"

for FILE in "$INDEX_DIR"/*.bam; do
    samtools index "$FILE"
done

echo "Indexing completed successfully. Files saved to $INDEX_DIR"
