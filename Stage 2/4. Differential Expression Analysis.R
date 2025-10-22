# DIFFERENTIAL EXPRESSION ANALYSIS - STAPHYLOCOCCUS AUREUS IN ACUTE VS CHRONIC PJIs

# Set working directory
setwd("C:/Users/HP/Dropbox/HackBio Internship  2025/Stage 2_RNA-Seq/Project_Microbes_S.aureus and PJIs/FeatureCounts and DEA")

# install/update packages and load libraries
install.packages("ggplot2")                          # For plotting
install.packages("ggfortify")                        # PCA plotting
install.packages("pheatmap")                         # For heatmaps
install.packages("BiocManager")                      # For installing Bioconductor packages
BiocManager::install("DESeq2")                       # For differential expression analysis
BiocManager::install("EnhancedVolcano")              # For volcano plots
library(ggplot2)                      
library(ggfortify)                    
library(pheatmap)                     
library(DESeq2)                       
library(EnhancedVolcano)
library(ggrepel)

# Set count file
count <- read.delim('counts.txt', header = T)
meta <- read.delim('metadata.tsv', header = T) 

# Preview
head(count)
head(meta)

# Keep to the important columns
raw_counts <- count[, meta$sample]
head(raw_counts)

# Add rownames
rownames(raw_counts) <- count$Geneid

head(raw_counts)

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = meta,
                              design = ~state)                  # Acute vs Chronic PJIs

# Preview
dds
dds$sample
dds$state

# Run DESeq2
dds <- DESeq(dds)

# Get results
res <- results(dds)

# Assess result
head(res)

# Filter for significant genes (adjusted p-value < 0.05)
sig_genes <- res[which(res$padj < 0.05), ]

# View result
print(sig_genes)

# PCA plot to visualize the separation of acute vs chronic samples
# First, check for constant columns
constant_columns <- apply(raw_counts, 2, function(x) length(unique(x)) <= 1)
constant_columns_indices <- which(constant_columns)
constant_columns_names <- colnames(raw_counts)[constant_columns_indices]

# View constant columns
print(constant_columns_names)
# Printed character(0)

# Check for zero variance rows
zero_variance_rows <- apply(raw_counts, 1, function(x) var(x) == 0)
zero_variance_row_indices <- which(zero_variance_rows)
zero_variance_row_names <- rownames(raw_counts)[zero_variance_row_indices]

# Print zero variance rows
print("Zero variance rows:")
print(zero_variance_row_names)

# Remove zero variance rows
raw_counts_filtered <- raw_counts[!zero_variance_rows, ]       
# Outcome: 2822 obs. of the 2942 obs. Filtered out 122 obs.

# Check for NA values
na_check <- any(is.na(raw_counts))
print(paste("Contains NA values:", na_check))
# Confirmed FALSE for the na_check.

# Now, create the PCA plot
pca_data <- prcomp(t(raw_counts_filtered), scale. = TRUE)
pca_df <- as.data.frame(pca_data$x)                                  # Gets PCA results               
pca_df$sample <- rownames(pca_df)                                    # Adds sample IDs
pca_df$state <- meta$state                                           # Adds metadata
explained_variance <- summary(pca_data)$importance[2, ] * 100        # Calculates percentage of variance explained by each principal component
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = state)) +   # Creates PCA plot
  geom_point(size = 3) +                                             # Points for each sample
  geom_text(aes(label = sample), vjust = 1.5, size = 3, check_overlap = TRUE) +  
  labs(
    title = "PCA of S. aureus Samples",
    x = paste("PC1 (", round(explained_variance[1], 2), "%)", sep = ""),
    y = paste("PC2 (", round(explained_variance[2], 2), "%)", sep = "")
  ) +
  theme_minimal() 
ggsave("pca_plot.png", plot = pca_plot)  # Saves PCA plot as a PNG file

# Create a volcano plot
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5, 5),
                ylim = c(0, -log10(1e-6)),
                title = 'Volcano Plot of DEGs',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)
ggsave("volcano_plot.png")                # Saves volcano plot as a PNG file

# Heatmap of significant DEGs
# Select significant genes for the heatmap
heatmap_data <- raw_counts[rownames(sig_genes), ]
pheatmap(heatmap_data, scale = "row", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         main = "Heatmap of Significant DEGs", show_rownames = TRUE)
ggsave("heatmap_plot.png")               # Save heatmap as a PNG file


# Identify the upregulated/downregulated genes
upregulated <- subset(res, padj < 0.05 & log2FoldChange > 2)
downregulated <- subset(res, padj < 0.05 & log2FoldChange < -2)

# List out the genes that are upregulated/downregulated
rownames(upregulated)
rownames(downregulated)

# Export the files
write.csv(upregulated, 'upregulated.csv')
write.csv(downregulated, 'downregulated.csv')
write.csv(raw_counts, 'raw_counts.csv')
write.csv(raw_counts_filtered, 'raw_counts_filtered.csv')

