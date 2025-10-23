# PROJECT 1: DIFFERENTIAL EXPRESSION ANALYSIS OF THE RESPONSE OF ARABIDOPSIS THALIANA VASCULATURE TISSUE TO UV TREATMENT
## INTRODUCTION
Environmental stresses affect plant tissues in various degrees, and plants' response to stress often have tissue-specific components. Berkowitz et al. (2021) previously compare the transcriptomes analysis of three leaf tissues (epidermis, mesophyll, vasculature) of _Arabidopsis thaliana_ in response to stress. This current project aims at determining the genes that respond to UV stress in the vasculature by: 
- performing differential expression analysis in the UV-treated vs water-treated vasculature (control);
- performing functional enrichment analysis and pathway mapping to identify major biological processes associated with the differentially expressed genes (DEGs).

## METHODS
### 1. Data acquisition
NCBI SRA database under project ID PRJNA668247

Replicate	Control	UV-C Treated
1	SRR12808527	SRR12808497
2	SRR12808528	SRR12808498
3	SRR12808529	SRR12808499

```bash
# Script for downloading
nano downloads.sh
```

```bash
#!/bin/bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/027/SRR12808527/SRR12808527.fastq.gz -o raw_data/SRR12808527_c1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/028/SRR12808528/SRR12808528.fastq.gz -o raw_data/SRR12808528_c2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/029/SRR12808529/SRR12808529.fastq.gz -o raw_data/SRR12808529_c3.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/097/SRR12808497/SRR12808497.fastq.gz -o raw_data/SRR12808497_t1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/098/SRR12808498/SRR12808498.fastq.gz -o raw_data/SRR12808498_t2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/099/SRR12808499/SRR12808499.fastq.gz -o raw_data/SRR12808499_t3.fastq.gz
```

```bash 
# Run script
bash downloads.sh
```

### 2. Preprocessing and quality control
#### 2a. Quality control with `fastqc`
```bash
# Write the script
nano quality_control.sh
```

```bash
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
    fastqc "$SAMPLE" -o "$OUTPUT_DIR" --quiet
done

echo "QC successfully done!"

# Aggregate QC reports with `multiqc`
echo "Aggregating QC reports..."

multiqc "$OUTPUT_DIR" -o "$MULTIQC_DIR"

echo "Reports saved to $MULTIQC_DIR"
```

```bash
# Run the script to perform QC
bash quality_control.sh
```

#### 2b. Trimming with `fastp`
```bash
# Write the script
nano trimming.sh
```

```bash
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
    --html "$TRIMMED_OUTPUT_DIR/${base}_fastp_report.html" --json "$TRIMMED_OUTPUT_DIR/${base}_fastp_report.json" \
    --overrepresented_threshold 0.01 \
    --cut_mean_quality 20 --cut_front --cut_tail \
    --qualified_quality_phred 20 --length_required 30 \
    --detect_adapter_for_pe 
done

echo "Results saved into $TRIMMED_OUTPUT_DIR." 
```

```bash
# Run the script to perform trimming
bash trimming.sh
```

#### 2c. Quality check on the trimmed data
```bash
# Write the script
nano trimmed_qc.sh
```

```bash
#!/bin/bash

# Set the variables
TRIMMED_OUTPUT_DIR="trimmed_reads"         # Directory containing the trimmed reads.
TRIMMED_QC_DIR="qc_trim"                   # Directory for trimmed reads QC reports.
TRIMMED_MULTIQC_DIR="multiqc_trim"         # Directory for the aggregated QC reports.

# Create output directories
mkdir -p "$TRIMMED_QC_DIR"
mkdir -p "$TRIMMED_MULTIQC_DIR"


echo "Performing quality check on the trimmed data..."

fastqc "$TRIMMED_OUTPUT_DIR"/*_trimmed.fastq -o "$TRIMMED_QC_DIR" --quiet
    
echo "Quality check successfully done!"

# Aggregate QC reports with `multiqc`
echo "Aggregating QC reports..."

multiqc "$TRIMMED_QC_DIR" -o "$TRIMMED_MULTIQC_DIR" --quiet

echo "Reports saved to $TRIMMED_MULTIQC_DIR"
```

```bash
# Run the script
bash trimmed_qc.sh
```

### 3. Mapping
```bash
# Create the genome directory
mkdir -p genome

# Change to the genome directory
cd genome/

# Download reference genome (TAIR10) from `Ensembl Plants`

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

# Unzip the zip file
gunzip Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz

#Rename the file for simplicity
mv Arabidopsis_thaliana.TAIR10.dna.toplevel.fa a_thaliana.fa

# create genome index directory
STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles a_thaliana.fa

# Return to the main working directory and progress with mapping
cd ../
```

```bash
# Script for mapping
nano mapping.sh
```

```bash
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
```

```bash
# Run script for mapping
bash mapping.sh
```

```bash
# Assess the bam output
samtools view mapped_reads/SRR12808497_t1_Aligned.sortedByCoord.out.bam | head -n 10       # Views the first 10 lines.
samtools view -c mapped_reads/SRR12808497_t1_Aligned.sortedByCoord.out.bam | head          # Counts total reads 
samtools flagstat mapped_reads/SRR12808497_t1_Aligned.sortedByCoord.out.bam | head         # Gets summary statistics for sample SRR20959676.

# Get summary statistics for all samples at a go.

for file in mapped_reads/*.bam; do
    echo "Processing $file"
    samtools flagstat "$file"
    echo "--------------------"
done
```

### 4. Counting the abundance of the transcriptome using `featureCounts`
```bash
# Download genome annotation
wget -nc -P genome https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz

# Unzip the gff3 file
gunzip ./genome/Arabidopsis_thaliana.TAIR10.62.gff3.gz

# Rename the gff3 file to a_thaliana.gff3 to be used for featureCounts
mv genome/Arabidopsis_thaliana.TAIR10.62.gff3 a_thaliana.gff3
```

```bash
# Write the script
nano featurecounts.sh
```

```bash
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
featureCounts -a $ANNOTATION_FILE -o $OUTPUT_COUNT_FILE -t gene -g gene_id -p -B -C -O "${BAM_FILES[@]}"

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
# -p: reads are paired-end
# -B: both reads of a pair map to the same feature
# -C: counts only reads that are uniquely mapped
# -O: counts reads overlapping multiple features only once - avoiding double-counting
```

```bash
# Run the script
bash featurecounts.sh
```

```bash
# Assess the featureCounts output 
cat counts/counts.txt | head -n 10      # First 10 lines.

# Download the counts.txt file into local device for downstream analyses in R Studio.
```


## REFERENCES
- Berkowitz, O., Xu, Y., Liew, L. C., Wang, Y., Zhu, Y., Hurgobin, B., Lewsey, M. G., and Whelan, J. (2021). RNA-seq analysis of laser microdissected _Arabidopsis thaliana_ leaf epidermis, mesophyll and vasculature defines tissue-specific transcriptional responses to multiple stress treatments. _The Plant Journal_. 104(3):938-955. DOI: 10.1111/tpj.15314.
