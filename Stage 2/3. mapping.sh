#! /bin/bash

# create directories for output
mkdir -p mapped

for infile in trimmed_reads/*.fastq ; do
    outfile=$(basename "$infile".fastq)
    STAR --genomeDir genome/genomeIndex --readFilesIn $infile --outFileNamePrefix mapped/$outfile --outSAMtype BAM SortedByCoordinate --outSAMattributes All
done