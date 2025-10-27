# DIFFERENTIAL EXPRESSION ANALYSIS - RESPONSE OF ARABIDOPSIS THALIANA VASCULATURE TISSUE TO STRESS

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


# FUNCTIONAL ENRICHMENT AND PATHWAY MAPPING

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



