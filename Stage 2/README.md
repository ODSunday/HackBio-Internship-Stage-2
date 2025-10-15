# PROJECT 2: TRANSCRIPTOMIC PROFILING OF _STAPHYLOCOCCUS AUREUS_ DURING ACUTE VS CHRONIC PHASES OF PERIPROSTHETIC JOINT INFECTION (PJI)
## INTRODUCTION
Periprosthetic joint infections (PJIs) are among the most devastating complications of orthopedic implants which increase morbidity, prolong hospital stays, and often require costly revision surgeries. Staphylococcus aureus, particularly methicillin-resistant strains (MRSA), is a leading cause of PJIs. 

A critical feature of _S. aureus_ is its ability to switch phenotypes between acute and chronic infection phases. For example, in **acute phase**, bacteria adopt an aggressive, planktonic growth mode, expressing virulence factors such as toxins, adhesins, and immune evasion genes, while, in **chronic phase**, bacteria adapt to a biofilm-like state, downregulating overt virulence and upregulating persistence pathways (stress response, metabolic rewiring, antibiotic tolerance). This adaptive flexibility makes chronic PJIs notoriously difficult to eradicate. Antibiotic regimens often fail, and host immune responses are blunted by biofilm shielding.

RNA sequencing uncovers the global transcriptional programmes underpinning the acute-to-chronic transition. Therefore, this project aims at capturing the gene expression profiles of _S. aureus_ isolates in the different clinical phases of PJI (acute and chronic) to: 
- Identify virulence, stress response, and metabolic genes that are significantly up- or downregulated;
- Perform functional enrichment analysis to highlight pathways associated with biofilm formation, immune evasion, and antibiotic resistance; and
- Explore small RNAs and regulatory elements potentially shaping the acute/chronic shift.

## METHODS
### 1. Data acquisition
The data for this project was obtained from SRA, using the SRA-Explorer:`https://sra-explorer.info/`. The following 7 samples in `PRJNA867318` were successfully downloaded into the `raw_data` directory. Each sample contains paired-end reads. 

|S/N|	Accession Number|	State|
|:---:|:-----:|:---------:|
|1|	SRR20959676| chronic periprosthetic joint infection|
|2|	SRR20959677| chronic periprosthetic joint infection|
|3|	SRR20959678| chronic periprosthetic joint infection|
|4|	SRR20959679| chronic periprosthetic joint infection|
|5|	SRR20959680| acute periprosthetic joint infection|
|6|	SRR20959681| acute periprosthetic joint infection|
|7|	SRR20959682| acute periprosthetic joint infection|

Script for downloading the datasets:
```bash
nano download.sh              # Writes the script containing the lines of command shown below:
```
```bash
#!/bin/bash

curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/077/SRR20959677/SRR20959677_1.fastq.gz -o SRR20959677_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/077/SRR20959677/SRR20959677_2.fastq.gz -o SRR20959677_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/076/SRR20959676/SRR20959676_1.fastq.gz -o SRR20959676_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/076/SRR20959676/SRR20959676_2.fastq.gz -o SRR20959676_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/082/SRR20959682/SRR20959682_1.fastq.gz -o SRR20959682_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/082/SRR20959682/SRR20959682_2.fastq.gz -o SRR20959682_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/081/SRR20959681/SRR20959681_1.fastq.gz -o SRR20959681_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/081/SRR20959681/SRR20959681_2.fastq.gz -o SRR20959681_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/080/SRR20959680/SRR20959680_1.fastq.gz -o SRR20959680_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/080/SRR20959680/SRR20959680_2.fastq.gz -o SRR20959680_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/079/SRR20959679/SRR20959679_1.fastq.gz -o SRR20959679_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/079/SRR20959679/SRR20959679_2.fastq.gz -o SRR20959679_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/078/SRR20959678/SRR20959678_1.fastq.gz -o SRR20959678_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/078/SRR20959678/SRR20959678_2.fastq.gz -o SRR20959678_2.fastq.gz
```
To download datasets:
```bash
bash download.sh
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
DATA_DIR="./raw_data"                        # Directory containing the raw data.
OUTPUT_DIR="./qc_reports"                    # Directory for the QC reports.
MULTIQC_DIR="./multiqc_reports"              # Directory for the aggregated QC reports.

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
DATA_DIR="./raw_data"                        # Directory containing the raw data.
TRIMMED_OUTPUT_DIR="./trimmed_reads"         # Directory for the trimmed reads.

# Create the output directory
mkdir -p "$TRIMMED_OUTPUT_DIR"

echo "Trimming reads with fastp..."

for FILE in "$DATA_DIR"/*_1.fastq.gz; do
    base=$(basename "$FILE" _1.fastq.gz)
    fastp -i "$DATA_DIR/${base}_1.fastq.gz" -I "$DATA_DIR/${base}_2.fastq.gz" \
    -o "$TRIMMED_OUTPUT_DIR/${base}_1.trimmed.fastq.gz" -O "$TRIMMED_OUTPUT_DIR/${base}_2.trimmed.fastq.gz" \
    --html "$TRIMMED_OUTPUT_DIR/${base}_fastp_report.html" --json "$TRIMMED_OUTPUT_DIR/${base}_fastp_report.json"
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
TRIMMED_OUTPUT_DIR="./trimmed_reads"         # Directory containing the trimmed reads.
TRIMMED_QC_DIR="./qc_trim"                   # Directory for trimmed reads QC reports.
TRIMMED_MULTIQC_DIR="./multiqc_trim"         # Directory for the aggregated QC reports.

# Create output directories
mkdir -p "$TRIMMED_QC_DIR"
mkdir -p "$TRIMMED_MULTIQC_DIR"

echo "Performing quality check on the trimmed data..."

fastqc "$TRIMMED_OUTPUT_DIR"/*.trimmed.fastq.gz -o "$TRIMMED_QC_DIR"
    
echo "Quality check successfully done!"

# Aggregate QC reports with `multiqc`
echo "Aggregating QC reports..."

multiqc "$TRIMMED_QC_DIR" -o "$TRIMMED_MULTIQC_DIR"

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

# Download reference genome from `Ensembl Bacteria`
wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-62/fasta/bacteria_41_collection/staphylococcus_aureus_gca_002895385/dna/Staphylococcus_aureus_gca_002895385.ASM289538v1_.dna.toplevel.fa.gz

# Unzip the zip file
gunzip Staphylococcus_aureus_gca_002895385.ASM289538v1_.dna.toplevel.fa.gz

#Rename the file for simplicity
mv Staphylococcus_aureus_gca_002895385.ASM289538v1_.dna.toplevel.fa s_aureus.fa
 
# create genome index directory
STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles s_aureus.fa

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
```
```bash
# Assess the bam output
samtools view mapped/SRR20959676_1.trimmed.fastq.fqAligned.sortedByCoord.out.bam | head
samtools view -c mapped/SRR20959676_1.trimmed.fastq.fqAligned.sortedByCoord.out.bam | head          
samtools flagstat mapped/SRR20959676_1.trimmed.fastq.fqAligned.sortedByCoord.out.bam | head

# Get summary statistics for all samples at a go.

for file in mapped/*.bam; do
    echo "Processing $file"
    samtools flagstat "$file"
    echo "--------------------"
done
```
