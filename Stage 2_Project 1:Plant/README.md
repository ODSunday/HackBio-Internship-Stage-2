# PROJECT 1: DIFFERENTIAL EXPRESSION ANALYSIS OF THE RESPONSE OF ARABIDOPSIS THALIANA VASCULATURE TISSUE TO UV TREATMENT
## INTRODUCTION
Environmental stresses affect plant tissues in various degrees, and plants' response to stress often have tissue-specific components. Berkowitz et al. (2021) previously compare the transcriptomes analysis of three leaf tissues (epidermis, mesophyll, vasculature) of _Arabidopsis thaliana_ in response to stress. This current project aims at determining the genes that respond to UV stress in the vasculature by: 
- performing differential expression analysis in the UV-treated versus water-treated vasculature (control);
- performing functional enrichment analysis and pathway mapping to identify major biological processes associated with the differentially expressed genes (DEGs).

## METHODS
### 1. Data acquisition
The datasets for this project (shown in the tables below) were from NCBI SRA database under project ID PRJNA668247, and obtained as follows:

###### Table 1: Samples
|Replicate|	Control|	UV-C Treated|
|:-------:|:------:|:--------------:|
|1|	SRR12808527	|SRR12808497|
|2|	SRR12808528	|SRR12808498|
|3|	SRR12808529	|SRR12808499|

###### Table 2: Renamed samples
|Replicate|	Control|	UV-C Treated|
|:-------:|:------:|:--------------:|
|1|	SRR12808527_c1	|SRR12808497_t1|
|2|	SRR12808528_c2	|SRR12808498_t2|
|3|	SRR12808529_c3	|SRR12808499_t3|

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
samtools view -c mapped_reads/SRR12808497_t1_Aligned.sortedByCoord.out.bam | head          # Counts total reads (26,900,918)
samtools flagstat mapped_reads/SRR12808497_t1_Aligned.sortedByCoord.out.bam | head         # Gets summary statistics for sample SRR12808497_t1.

# Get summary statistics for all samples at a go.

for file in mapped_reads/*.bam; do
    echo "Processing $file"
    samtools flagstat "$file"
    echo "--------------------"
done
```

### 4. Visualizing the Transcriptome Map

```bash
# Create a directory for IGV viewing
mkdir -p igv

# Copy the mapped reads (output BAM files) into the igv directory for IGV viewing
cp mapped_reads/*.bam igv/
```

#### 4a. Indexing
```bash
# Script for indexing
nano igv_indexing.sh
```

```bash
#!/bin/bash

INDEX_DIR="./igv"

echo "Indexing in progress"
echo "--------------------"

for FILE in "$INDEX_DIR"/*.bam; do
    samtools index "$FILE"
done

echo "Indexing completed successfully. Files saved to $INDEX_DIR"
```

```bash
# Run script for indexing
bash igv_indexing.sh

# Download the igv directory to the local computer for visualizing in IGV 
```

#### 4b. Visualize in IGV

`https://igv.org/app/`

### 5. Counting the abundance of the transcriptome using `featureCounts`
```bash
# Download genome annotation
wget -nc -P genome https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.62.gff3.gz

# Unzip the gff3 file
gunzip ./genome/Arabidopsis_thaliana.TAIR10.62.gff3.gz

# Rename the gff3 file to a_thaliana.gff3 to be used for featureCounts
mv ./genome/Arabidopsis_thaliana.TAIR10.62.gff3 ./genome/a_thaliana.gff3
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

### 6. Differential expression analysis, functional enrichment, and pathway mapping: `R Studio`

```r
## DIFFERENTIAL EXPRESSION ANALYSIS

# Setting working directory
setwd("C:/Users/HP/Dropbox/HackBio Internship  2025/Stage 2_RNA-Seq/Project_Plant_Arabidopsis thaliana")

# installing/updating packages and loading libraries
install.packages("ggplot2")                          # For plotting
install.packages("ggfortify")                        # PCA plotting
install.packages("pheatmap")                         # For heatmaps
install.packages("dplyr")                            # For data manipulation
install.packages("BiocManager")                      # For installing Bioconductor packages
BiocManager::install("DESeq2")                       # For differential expression analysis
BiocManager::install("EnhancedVolcano")              # For volcano plots
BiocManager::install("clusterProfiler")              # For GO/KEGG enrichment analysis
BiocManager::install("org.At.tair.db")               # For annotating Arabidopsis thaliana
BiocManager::install("enrichplot")                   # For visualizing enrichment results
library(ggplot2)                      
library(ggfortify)                    
library(pheatmap)
library(dplyr)
library(DESeq2)                       
library(EnhancedVolcano)
library(ggrepel)
library(clusterProfiler)
library(org.At.tair.db)    
library(enrichplot)                   

# Setting count file
count <- read.delim('counts.txt', header = T)
meta <- read.delim('metadata.tsv', header = T) 

# Preview
head(count)
head(meta)

# Keepping to the important columns
raw_counts <- count[, meta$sample]
head(raw_counts)

# Adding rownames
rownames(raw_counts) <- count$Geneid

head(raw_counts)

# Creating DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = meta,
                              design = ~ condition)       # uv-treated vs control                

# Removing rows with low counts
dds <- dds[rowSums(counts(dds)) > 1, ]

# Preview
dds
dds$sample
dds$condition

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Assess result
head(res)

# Summary of results
summary(res)

# Using dplyr order results by adjusted p-value and convert to data frame
resOrdered <- as.data.frame(res) %>%
  arrange(padj)  

# Save the results to a CSV file
write.csv(resOrdered, file = "DE_results.csv")

# List of top 100 differentially expressed genes (DEGs)
top100 <- resOrdered %>%
  filter(!is.na(padj)) %>%
  top_n(-100, padj)                                # Gets the top 100 based on adjusted p-value
write.csv(top100, file = "top100_DEGs.csv")

# Creating MA plot
plotMA(res, ylim = c(-5, 5), main = "MA Plot")
ggsave("MA_plot.png")

# Creating PCA to visualize sample clustering
rld <- rlog(dds)
pcaData <- plotPCA(rld, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Generating the PCA plot
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  geom_text(aes(label = rownames(pcaData)), vjust = -1, hjust = 1, size = 3) +  # To add sample names
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of UV-C Treated vs Control") +
  theme_minimal()
ggsave("PCA_plot1.png")      # Saves the PCA plot              

# Creating a heatmap for the top 20 DEGs
top20 <- resOrdered %>%
  filter(!is.na(padj)) %>%
  top_n(-20, padj)  # Get top 20 based on adjusted p-value
mat <- assay(rld)[rownames(top20), ]
mat <- mat - rowMeans(mat)                    # Centres the data
pheatmap(mat, cluster_rows = TRUE, cluster_cols = TRUE, 
         show_rownames = TRUE, 
         main = "Heatmap of Top 20 DEGs")
ggsave("heatmap.png")                       # Saves the heatmap

# Creating an Enhanced Volcano plot
p <- EnhancedVolcano(resOrdered,
                     lab = rownames(resOrdered),
                     x = 'log2FoldChange',
                     y = 'padj',
                     xlim = c(-6, 6),                                       # Sets the x-axis limits
                     ylim = c(0, 10),                                       # Sets the y-axis limits
                     pCutoff = 0.05,                                        # Sets p-value cutoff
                     FCcutoff = 1,                                          # Sets fold change cutoff
                     pointSize = 3.0,                                       # Sets size of points
                     labSize = 3.0,                                         # Sets size of labels 
                     caption = 'Red points repreent significant genes',     # Adds caption
                     col = c('grey30', 'forestgreen', 'royalblue', 'red3'), # Adds colours for points
                     legendLabels = c('Not significant', 'Log2 FC > 1', 'Log2 FC < -1', 'Adjusted p < 0.05'),  # Legend labels
                     legendPosition = "top",                                # Moves legend to the top
                     legendLabSize = 10,                                    # Increases legend label size for better visibility
                     legendIconSize = 4,                                    # Adjusts the size of the legend icons
                     drawConnectors = TRUE,                                 # Draws connectors for better visibility
                     arrowheads = TRUE)                                     # Adds arrowheads to connectors
p <- p + 
  ggtitle("Volcano Plot", "Differential Expression Analysis") +             # Adds title and subtitle
  theme(plot.title = element_text(hjust = 0.5),                             # Centres the title
        plot.subtitle = element_text(hjust = 0.5))                          # Centres the subtitle
ggsave("enhanced_volcano_plot.png", plot = p, width = 12, height = 8, dpi = 300)  # Saves the plot

## FUNCTIONAL ENRICHMENT AND PATHWAY MAPPING

# Getting the top 5 upregulated and downregulated genes
upregulated <- rownames(resOrdered[resOrdered$log2FoldChange > 1 & !is.na(resOrdered$padj), ])
downregulated <- rownames(resOrdered[resOrdered$log2FoldChange < -1 & !is.na(resOrdered$padj), ])

## A: Gene Ontology (GO) Enrichment Analysis
# Performing enrichment analysis for upregulated genes across all GO categories
up_enrich_BP <- enrichGO(gene = upregulated,
                      OrgDb = org.At.tair.db,
                      keyType = "TAIR",
                      ont = "BP",                            # Biological Process
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

up_enrich_CC <- enrichGO(gene = upregulated,
                         OrgDb = org.At.tair.db,
                         keyType = "TAIR",
                         ont = "CC",                         # Cellular Component
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05)

up_enrich_MF <- enrichGO(gene = upregulated,
                         OrgDb = org.At.tair.db,
                         keyType = "TAIR",
                         ont = "MF",                         # Molecular Function
                         pAdjustMethod = "BH",
                         qvalueCutoff = 0.05)

# Performing enrichment analysis for downregulated genes across all GO categories
down_enrich_BP <- enrichGO(gene = downregulated,
                        OrgDb = org.At.tair.db,
                        keyType = "TAIR",
                        ont = "BP",                          # Biological Process
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

down_enrich_CC <- enrichGO(gene = downregulated,
                           OrgDb = org.At.tair.db,
                           keyType = "TAIR",
                           ont = "CC",                       # Cellular Component
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05)

down_enrich_MF <- enrichGO(gene = downregulated,
                           OrgDb = org.At.tair.db,
                           keyType = "TAIR",
                           ont = "MF",                       # Molecular Function
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05)

# Viewing the results for upregulated genes
head(up_enrich_BP)
head(up_enrich_CC)
head(up_enrich_MF)

# Viewing the results for downregulated genes
head(down_enrich_BP)
head(down_enrich_CC)
head(down_enrich_MF)

# Visualizing enrichment results for upregulated genes
dot_plot <- dotplot(up_enrich_BP, showCategory = 10) +       # Creates the dot plot
  ggtitle("Top 10 Enriched GO Terms (Upregulated Genes - BP)") +           
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   # Rotates x-axis labels
        plot.margin = margin(10, 10, 10, 10))                # Adjusts margins
print(dot_plot)                                              # Prints the plot
ggsave("dotplot_upregulated_BP.png", plot = dot_plot, width = 10, height = 6, dpi = 300) # Saves the plot

dot_plot <- dotplot(up_enrich_CC, showCategory = 10) +       
  ggtitle("Top 10 Enriched GO Terms (Upregulated Genes - CC)") +           
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   
        plot.margin = margin(10, 10, 10, 10))                
print(dot_plot)                                              
ggsave("dotplot_upregulated_CC.png", plot = dot_plot, width = 10, height = 6, dpi = 300) 

dot_plot <- dotplot(up_enrich_MF, showCategory = 10) +       
  ggtitle("Top 10 Enriched GO Terms (Upregulated Genes - MF)") +           
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   
        plot.margin = margin(10, 10, 10, 10))                
print(dot_plot)                                              
ggsave("dotplot_upregulated_MF.png", plot = dot_plot, width = 10, height = 6, dpi = 300) 

# Visualizing enrichment results for downregulated genes
dot_plot <- dotplot(down_enrich_BP, showCategory = 10) +     # Creates the dot plot
  ggtitle("Top 10 Enriched GO Terms (Downregulated Genes - BP)") +           
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   # Rotates x-axis labels
        plot.margin = margin(10, 10, 10, 10))                # Adjusts margins
print(dot_plot)                                              # Prints the plot
ggsave("dotplot_downregulated_BP.png", plot = dot_plot, width = 10, height = 6, dpi = 300) # Saves the plot

dot_plot <- dotplot(down_enrich_CC, showCategory = 10) +     
  ggtitle("Top 10 Enriched GO Terms (Downregulated Genes - CC)") +           
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   
        plot.margin = margin(10, 10, 10, 10))                
print(dot_plot)                                              
ggsave("dotplot_downregulated_CC.png", plot = dot_plot, width = 10, height = 6, dpi = 300) 

dot_plot <- dotplot(down_enrich_MF, showCategory = 10) +     
  ggtitle("Top 10 Enriched GO Terms (Downregulated Genes - MF)") +           
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   
        plot.margin = margin(10, 10, 10, 10))                
print(dot_plot)                                              
ggsave("dotplot_downregulated_MF.png", plot = dot_plot, width = 10, height = 6, dpi = 300) 

# Getting the top 5 enriched pathways for upregulated genes
top5_up_BP <- head(up_enrich_BP, n = 5)
top5_up_CC <- head(up_enrich_CC, n = 5)
top5_up_MF <- head(up_enrich_MF, n = 5)

# Printing the results of the top 5 enriched pathways for upregulated genes
print("Top 5 Enriched Biological Processes (Upregulated Genes):")
print(top5_up_BP)

print("Top 5 Enriched Cellular Components (Upregulated Genes):")
print(top5_up_CC)

print("Top 5 Enriched Molecular Functions (Upregulated Genes):")
print(top5_up_MF)

# Getting the top 5 enriched pathways for downregulated genes
top5_down_BP <- head(down_enrich_BP, n = 5)
top5_down_CC <- head(down_enrich_CC, n = 5)
top5_down_MF <- head(down_enrich_MF, n = 5)

# Printing the results of the top 5 enriched pathways for downregulated genes
print("Top 5 Enriched Biological Processes (Downregulated Genes):")
print(top5_down_BP)

print("Top 5 Enriched Cellular Components (Downregulated Genes):")
print(top5_down_CC)

print("Top 5 Enriched Molecular Functions (Downregulated Genes):")
print(top5_down_MF)

# Summarizing the results of GO enrichment analysis
summarize_GO <- function(enrichment_results, category) {
  summary <- data.frame(
    Term = enrichment_results$Description,
    GeneCount = enrichment_results$Count,
    PValue = enrichment_results$pvalue,
    QValue = enrichment_results$qvalue,
    Category = category
  )
  return(summary)
}

# Summarizing the top 5 enriched pathways for upregulated genes
summary_up_BP <- summarize_GO(top5_up_BP, "Biological Process")
summary_up_CC <- summarize_GO(top5_up_CC, "Cellular Component")
summary_up_MF <- summarize_GO(top5_up_MF, "Molecular Function")

# Summarizing the top 5 enriched pathways for downregulated genes
summary_down_BP <- summarize_GO(top5_down_BP, "Biological Process")
summary_down_CC <- summarize_GO(top5_down_CC, "Cellular Component")
summary_down_MF <- summarize_GO(top5_down_MF, "Molecular Function")

# Combining all summary results into a single data frame
summary_go_results <- rbind(summary_up_BP, summary_up_CC, summary_up_MF,
                         summary_down_BP, summary_down_CC, summary_down_MF)

# Printing the summary results
print(summary_go_results)

write.csv(summary_go_results, "go_enrichment_summary_results.csv", row.names = FALSE)

## B: KEGG Enrichment Analysis
# Performing KEGG enrichment analysis for upregulated genes
up_kegg <- enrichKEGG(gene = upregulated,
                      organism = 'ath',               # For Arabidopsis thaliana
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

# Performing KEGG enrichment analysis for downregulated genes
down_kegg <- enrichKEGG(gene = downregulated,
                        organism = 'ath',  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

# Visualizing the KEGG results (top 10) for upregulated genes
dot_plot <- dotplot(up_kegg, showCategory = 10) + ggtitle("KEGG Enrichment for Upregulated Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   
      plot.margin = margin(10, 10, 10, 10))                
print(dot_plot)
ggsave("upregulated_kegg_dotplot.png", plot = dot_plot, width = 10, height = 6, dpi = 300)

# Visualizing the KEGG results (top 10) for downregulated genes
dot_plot <- dotplot(down_kegg, showCategory = 10) + ggtitle("KEGG Enrichment for Downregulated Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),   
        plot.margin = margin(10, 10, 10, 10))                
print(dot_plot)
ggsave("downregulated_kegg_dotplot.png", plot = dot_plot, width = 10, height = 6, dpi = 300)

# Printing the results for upregulated genes
print("KEGG Enrichment Results for Upregulated Genes:")
head(up_kegg)

# Printing the results for downregulated genes
print("KEGG Enrichment Results for Downregulated Genes:")
head(down_kegg)

# Printing the top 5 KEGG pathways for upregulated genes
top5_up_kegg <- head(up_kegg, n = 5)
print("Top 5 KEGG Pathways for Upregulated Genes:")
print(top5_up_kegg)

# Printing the top 5 KEGG pathways for downregulated genes
top5_down_kegg <- head(down_kegg, n = 5)
print("Top 5 KEGG Pathways for Downregulated Genes:")
print(top5_down_kegg)

# Summarizing the results for KEGG enrichment analysis
summarize_kegg <- function(kegg_results, category) {
  summary <- data.frame(
    Pathway = kegg_results$Description,
    GeneCount = kegg_results$Count,
    PValue = kegg_results$pvalue,
    QValue = kegg_results$qvalue,
    Category = category
  )
  return(summary)
}

# Summarizing the top 5 KEGG pathways for upregulated genes
summary_up_kegg <- summarize_kegg(top5_up_kegg, "Upregulated")

# Summarizing the top 5 KEGG pathways for downregulated genes
summary_down_kegg <- summarize_kegg(top5_down_kegg, "Downregulated")

# Combining all summary results into a single data frame
summary_kegg_results <- rbind(summary_up_kegg, summary_down_kegg)

# Printing the summary results
print(summary_kegg_results)

write.csv(summary_kegg_results, "kegg_enrichment_summary_results.csv", row.names = FALSE)
```

## RESULTS AND DISCUSSION
### Principal component analysis (PCA) shows a clear separation between treatment and control
PCA revealed a clear separation between UV-C treated and control vascular tissue samples of _Arabidopsis thaliana_ (**Figure 1**). The first principal component (PC1) explained 94% of the total variance, while the second principal component (PC2) explained 2%, indicating that the main source of variation was due to UV-C treatment. The clear separation between the UV-C treated samples ('t1–'t3) and the control samples ('c1–'c3) along PC1 indicates distinct transcriptional profiles between the two conditions. Hence, the PCA result confirms that UV-C exposure induces significant transcriptional reprogramming in the vascular tissues, likely involving stress-responsive, DNA repair, and antioxidant pathways (Berkowitz et al. 2021), thereby providing a solid basis for subsequent differential expression and pathway enrichment analyses to identify specific genes and molecular pathways responsive to UV-C stress.

<img width="700" height="600" alt="PCA_plot1" src="https://github.com/user-attachments/assets/43cc4392-bccc-4ff1-8898-fc6eb1322c59" />

###### Figure 1: Principal component analysis showing the clustering of UV-C treated samples compared to control

### UV stress induced significant transcriptional differences in differentially expressed genes (DEGs) 
Apart from filtering the top 100 DEGs, a heatmap of the top 20 DEGs was generated to visualize the transcriptional differences between UV-C treated and control vascular tissues of _A. thaliana_ (**Figure 2**). The heatmap shows two distinct clusters corresponding to the experimental groups. The UV-C treated samples ('t1–'t3) grouped together and displayed predominantly higher expression levels (red) compared to the control samples ('c1–'c3) which formed a separate cluster with lower expression values (blue). 

This clear clustering pattern is consistent with the PCA results, confirming the strong differential gene expression response induced by UV-C exposure. The upregulated genes in the UV-treated samples may be involved in UV-protective mechanisms such as DNA repair, oxidative stress defense, and secondary metabolite biosynthesis—processes commonly activated under UV stress conditions. Conversely, the downregulated genes in the control samples likely represent pathways suppressed by UV exposure, possibly linked to normal growth and metabolic activities that are not prioritised under stress.

Hence, the heatmap demonstrates that UV-C stress triggers a coordinated transcriptional reprogramming in the vascular tissues, resulting in the distinct upregulation of stress-responsive genes and the downregulation of non-essential pathways. These patterns provide a molecular basis for understanding how the vascular system of _A. thaliana_ adapts to UV-induced damage and maintains physiological stability under abiotic stress conditions.

Furthermore, a volcano plot was generated to visualize the overall distribution of DEGs between UV-C treated and control vascular tissues of Arabidopsis thaliana (**Figure 3**). Each point represents an individual gene, plotted according to its log₂ fold change (x-axis) and statistical significance (–log₁₀ p-value, y-axis).

The plot shows a large number of significantly regulated genes (adjusted p < 0.05), highlighted in red. Genes with positive log₂ fold change values are upregulated under UV-C treatment, while those with negative values are downregulated relative to the control. The clear asymmetry and wide spread of significant points along both directions indicate that UV-C exposure triggers substantial transcriptional reprogramming in the vascular tissues, involving both activation and repression of numerous genes.

Notably, several genes exhibited more than 4-fold changes in expression, underscoring strong UV-induced responses. These genes likely participate in photoprotective and stress-response pathways, including DNA repair, oxidative stress mitigation, and secondary metabolism. The volcano plot thus provides a comprehensive overview of the extent and magnitude of transcriptional alterations driven by UV-C treatment, reinforcing the PCA and heatmap findings that demonstrated distinct expression patterns between the treated and control groups.

<img width="655" height="768" alt="heatmap" src="https://github.com/user-attachments/assets/f9b0c32b-fadc-452e-b5d8-f610aa96730d" />

###### Figure 2: Heatmap of top 20 DEGS

<img width="3600" height="2400" alt="enhanced_volcano_plot" src="https://github.com/user-attachments/assets/e5af33c0-53f5-4aa4-86c7-a1ecd3e21bca" />

###### Figure 3: Volcano plot of the differential exprssion analysis showing upregulated and downregulated DEGs 

### Functional enrichment and pathway analyses of the top upregulated and downregulated genes
Given the extensive transcriptional alterations observed in the UV-C treated vascular tissues, functional enrichment analysis and pathway mapping were subsequently performed to identify the major biological processes and molecular pathways affected by the treatment. The aim was to determine whether the significantly upregulated and downregulated genes are functionally associated with specific stress-response mechanisms to provide more insights into the molecular basis of UV-induced adaptations in _A. thaliana_ vascular tissues.

First, Gene Ontology (GO) enrichment analysis was performed to gain functional insights into the biological roles of genes significantly affected by UV-C exposure in _A. thaliana_ vascular tissues. The top enriched GO terms were observed for both upregulated and downregulated genes across the three major GO categories: Biological Process (BP), Cellular Component (CC), and Molecular Function (MF) (**Table 3**; **Figures 4** and **5**).

**Top GO enriched terms and pathways for upregulated genes**

The upregulated genes were predominantly enriched in stress-related biological processes, particularly those associated with cellular response to decreased oxygen levels, response to hypoxia, and response to oxygen levels. This pattern suggests that UV-C treatment triggers hypoxia-like stress responses in vascular tissues, possibly due to oxidative damage or disrupted cellular respiration. The enrichment of terms related to the nucleolus, preribosome, and cytoplasmic stress granule in the Cellular Component category indicates active ribosome biogenesis and translational regulation under stress conditions.

At the molecular level, upregulated genes showed strong enrichment for calmodulin binding, glutathione transferase activity, and transmembrane receptor protein kinase activity, implying the involvement of calcium-mediated signaling, detoxification pathways, and receptor-mediated stress signaling. Collectively, these results highlight the activation of defensive and adaptive pathways aimed at maintaining redox balance, regulating protein synthesis, and enhancing cellular resilience under UV-C stress.

**GO enriched terms and pathways for downregulated genes**

Conversely, the downregulated genes were enriched in biological processes linked to plastid and chloroplast organization, RNA modification, and cell cycle regulation. These findings suggest that UV-C exposure suppresses photosynthetic and growth-related activities, likely as an energy-conservation strategy during stress. The corresponding Cellular Component terms (chloroplast thylakoid, plastid nucleoid, and chloroplast organization) further reinforce the downregulation of photosynthetic machinery.

In terms of molecular function, the enrichment of kinase regulator activity and cyclin-dependent protein serine/threonine kinase regulator activity points to reduced cell cycle progression and signaling activity. This suppression may reflect a cellular shift from growth-oriented to stress-protective metabolism following UV-C exposure.

Summarily, the GO enrichment results reveal that UV-C treatment induces a dual response in vascular tissues: activation of stress-responsive, redox-regulatory, and signaling pathways, and repression of photosynthetic and cell cycle processes. This coordinated transcriptional reprogramming underscores the plant’s adaptive strategy to mitigate UV-induced damage while maintaining essential metabolic balance.

**A** <img width="1000" height="600" alt="dotplot_upregulated_BP" src="https://github.com/user-attachments/assets/1742f88f-44cb-4183-9e96-4aa4f201710e" />

**B** <img width="1000" height="600" alt="dotplot_upregulated_CC" src="https://github.com/user-attachments/assets/8e1ed0f2-639f-4f06-b1bb-0646309e6171" />

**C** <img width="1000" height="600" alt="dotplot_upregulated_MF" src="https://github.com/user-attachments/assets/4dde5cf8-a730-4583-8993-6a6c3bdc68b6" />

###### Figures 4: Top 10 enriched GO terms for the upregulated genes. (A): Biological process (BP). (B): Cellular components (CC). (C): Molecular function (MF).

**A** <img width="1000" height="600" alt="dotplot_downregulated_BP" src="https://github.com/user-attachments/assets/227527d4-43bb-4fc8-8358-99b9fd0b7181" />

**B** <img width="1000" height="600" alt="dotplot_downregulated_CC" src="https://github.com/user-attachments/assets/36e6c6ab-701e-4ef0-b801-3ae7028f9572" />

**C** <img width="1000" height="600" alt="dotplot_downregulated_MF" src="https://github.com/user-attachments/assets/992d86e6-4938-46df-b31b-3b95ff6c95b1" />

###### Figures 5: Top 10 enriched GO terms for the downregulated genes. (A): Biological process (BP). (B): Cellular components (CC). (C): Molecular function (MF).

###### Table 3: Top 5 GO enriched pathways for the upregulated and downregulated genes

|Term|	Gene Count|	Category|
|:--:|:----------:|:-------:|
**Upregulated DEGs**				
|cellular response to decreased oxygen levels|	153|	Biological Process|
|cellular response to oxygen levels	|153	|Biological Process|
|cellular response to hypoxia|	152|	Biological Process|
|response to oxygen levels|	163	|Biological Process|
|response to hypoxia|	161|	Biological Process|
|preribosome|	38|	Cellular Component|
|nucleolus	|130	|Cellular Component|
|90S preribosome	|15	|Cellular Component|
|preribosome, large subunit precursor|	12|	Cellular Component|
|cytoplasmic stress granule|	49|	Cellular Component|
|calmodulin binding	|74	|Molecular Function|
|glutathione transferase activity	|27	|Molecular Function|
|ADP binding	|57|	Molecular Function|
|transmembrane signaling receptor activity|	37	|Molecular Function|
|transmembrane receptor protein kinase activity|	28	|Molecular Function|
**Downregulated DEGs**				
|plastid organization	|115|	Biological Process|
|chloroplast organization|	92|	Biological Process|
|RNA modification	|114	|Biological Process|
|mitotic cell cycle|	75	|Biological Process|
|cell cycle phase transition|	39|	Biological Process|
|chloroplast nucleoid|	31	|Cellular Component|
|plastid nucleoid	|32|	Cellular Component|
|chloroplast thylakoid|	127	|Cellular Component|
|nucleoid	|33	|Cellular Component|
|plastid thylakoid|	127|	Cellular Component|
|cyclin-dependent protein serine/threonine kinase regulator activity|	25|	Molecular Function|
|kinase regulator activity	|37|	Molecular Function|
|protein kinase regulator activity|	35|	Molecular Function|

#### Kyoto Encyclopedia of Genes and Genomes (KEGG) enrichment analysis for the upregulated and downregulated genes
To further understand the biological pathways affected by UV-C exposure, KEGG pathway enrichment analysis was performed for the significantly upregulated and downregulated genes in the vascular tissues of _A. thaliana_ (**Table 4**; **Figure 6**).

**Top KEGG enriched pathways for upregulated genes**

The upregulated genes were significantly enriched in pathways associated with plant–pathogen interaction, MAPK signaling, and phenylpropanoid biosynthesis. Enrichment of the plant–pathogen interaction and MAPK signaling pathways suggests activation of defense-related signaling cascades that mediate stress perception and response. The phenylpropanoid biosynthesis pathway, known for producing antioxidant and UV-absorbing secondary metabolites (e.g., flavonoids and lignin precursors), indicates a protective response aimed at minimizing UV-induced oxidative damage.

Additionally, enrichment in glutathione metabolism underscores the role of redox homeostasis and detoxification mechanisms in combating reactive oxygen species generated under UV-C stress. The ribosome biogenesis in eukaryotes pathway suggests enhanced protein synthesis capacity to support repair and defense processes. Together, these results imply that UV-C exposure activates important protective mechanisms involving stress signaling, antioxidant defense, and metabolic reprogramming.

**Top KEGG enriched pathways for downregulated genes**

The downregulated genes were primarily enriched in biosynthesis of cofactors, glyoxylate and dicarboxylate metabolism, and glycine, serine, and threonine metabolism, indicating suppression of core metabolic and biosynthetic activities. The reduced expression of genes in the plant hormone signal transduction pathway suggests altered hormonal balance and signaling, potentially contributing to growth inhibition and stress adaptation. Similarly, the downregulation of folate biosynthesis may reflect a redirection of metabolic resources away from growth and toward stress defense.

The KEGG enrichment analysis complements the GO results, demonstrating that UV-C treatment reprogrammes vascular tissue metabolism by upregulating stress-related and defense-related pathways while downregulating pathways involved in growth, photosynthesis, and biosynthesis. This coordinated shift underscores the plant’s adaptive strategy to withstand UV-induced oxidative and metabolic stress.

**A**
<img width="1000" height="600" alt="upregulated_kegg_dotplot" src="https://github.com/user-attachments/assets/25f1f496-bae2-416d-8f1e-ecf341dad765" />

**B**
<img width="1000" height="600" alt="downregulated_kegg_dotplot" src="https://github.com/user-attachments/assets/6346225c-1248-4b02-be29-18778b3aba36" />

###### Figure 6: Top 10 enriched KEGG pathways for the DEGs. (A) Upregulated DEGs. (B) Downregulated DEGs

###### Table 4: Top 5 KEGG enriched pathways for the upregulated and downregulated genes
|Pathway	|Gene Count	|Category|
|:---------:|:---------:|:------:|
|Plant-pathogen interaction	|93|	Upregulated|
|MAPK signaling pathway - plant|	46	|Upregulated|
|Phenylpropanoid biosynthesis|	42|	Upregulated|
|Glutathione metabolism	|33	|Upregulated|
|Ribosome biogenesis in eukaryotes|	32	|Upregulated|
|Biosynthesis of cofactors	|60	|Downregulated|
|Glyoxylate and dicarboxylate metabolism	|27	|Downregulated|
|Glycine, serine and threonine metabolism	|24	|Downregulated|
|Plant hormone signal transduction	|91	|Downregulated|
|Folate biosynthesis|	14	|Downregulated|

## CONCLUSION
UV exposure induces broad and coordinated changes in gene expression in the vascular tissues of _A. thaliana_. The UV-C treatment upregulated stress-related pathways such as hypoxia response, redox regulation, and signaling activity, while downregulating genes associated with chloroplast organization and cell cycle progression. Similarly, UV exposure activates defense and detoxification mechanisms, particularly plant–pathogen interaction, MAPK signaling, phenylpropanoid biosynthesis, and glutathione metabolism, while suppressing metabolic and hormone signaling pathways. Hence, vascular tissues of _A. thaliana_ respond to UV-C stress through an integrated network of defense signaling, antioxidant activity, and metabolic reprogramming targeted at maintaining cellular homeostasis and enhancing stress resilience.

## REFERENCES
- Berkowitz, O., Xu, Y., Liew, L. C., Wang, Y., Zhu, Y., Hurgobin, B., Lewsey, M. G., and Whelan, J. (2021). RNA-seq analysis of laser microdissected _Arabidopsis thaliana_ leaf epidermis, mesophyll and vasculature defines tissue-specific transcriptional responses to multiple stress treatments. _The Plant Journal_. 104(3):938-955. DOI: 10.1111/tpj.15314.
