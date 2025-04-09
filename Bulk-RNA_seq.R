###############################
### 1. Load Required Libraries
###############################
# (Uncomment installation lines if packages are not installed)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "org.Hs.eg.db", "DESeq2", "apeglm", "pheatmap", "vsn"))
# install.packages(c("dplyr", "RColorBrewer", "readr", "openxlsx", "factoextra"))

library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(vsn)
library(RColorBrewer)
library(readr)
library(openxlsx)
library(factoextra)
library(tibble)
library(openxlsx)

#####################################
### 2. Data Loading and Preprocessing
#####################################
setwd("C:\\Users\\hp\\Desktop\\R studio")
kimdata <- read.csv("C:\\Users\\hp\\Desktop\\R studio\\update code with STAT1 filter\\Thomas.rawcounts.csv")

# Define treatment labels and assign as column names
treatment <- c("X", "48h LPS", "48h C", "72h LPS", "48h C Rep1",
               "24h C", "72h C", "24h C Rep1",
               "24h LPS", "72h C Rep1", "48h LPS Rep1",
               "24h C Rep2", "24h LPS Rep1", "24h LPS Rep2",
               "72h LPS Rep1", "48h LPS Rep2", "48h C Rep2")
colnames(kimdata) <- treatment

# Reorder columns for consistency
kimdata <- kimdata[, c("X", "24h C", "24h C Rep1", "24h C Rep2", 
                       "24h LPS", "24h LPS Rep1", "24h LPS Rep2",
                       "48h C", "48h C Rep1", "48h C Rep2", 
                       "48h LPS", "48h LPS Rep1", "48h LPS Rep2",
                       "72h C", "72h C Rep1", "72h LPS", "72h LPS Rep1")]

# Convert Ensembl IDs to gene symbols
geneID <- select(org.Hs.eg.db, keys = kimdata$X, keytype = "ENSEMBL", columns = "SYMBOL")
geneID_unique <- geneID[!duplicated(geneID[c("ENSEMBL")]), ]
geneID_clean <- na.omit(geneID_unique)
Symbols_old <- geneID_unique$SYMBOL  # May contain NA
Symbols <- geneID_clean$SYMBOL         # Clean list without NA

# Create a counts matrix (remove original Ensembl column)
counts_data <- subset(kimdata, select = -c(X))
# Use coalesce to fill in gene names and ensure uniqueness
Rows <- dplyr::coalesce(Symbols_old, rownames(counts_data))
Unique.Rows <- make.unique(Rows, sep = " ")
rownames(counts_data) <- Unique.Rows
counts_data <- counts_data[rownames(counts_data) %in% Symbols, ]

# Define experimental conditions and time points
condition <- c("untreated", "untreated", "untreated",
               "treated", "treated", "treated",
               "untreated", "untreated", "untreated",
               "treated", "treated", "treated", 
               "untreated", "untreated", "treated", "treated")
time <- c("24h", "24h", "24h",
          "24h", "24h", "24h", 
          "48h", "48h", "48h",
          "48h", "48h", "48h", 
          "72h", "72h", "72h", "72h")
samples <- data.frame(condition = condition, time = time, row.names = colnames(counts_data))

##########################################
### 3. Filtering and Differential Analysis
##########################################
# Remove genes with all-zero counts
counts_filtered <- counts_data %>%
  filter(rowSums(across(where(is.numeric))) != 0)

# Create DESeq2 dataset and run analysis
dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = samples,
                              design = ~ condition + time)
dds <- DESeq(dds)
# Keep genes with a total count of at least 10
dds <- dds[rowSums(counts(dds)) >= 10, ]
ntd <- normTransform(dds)

# Export filtered counts for record
write.xlsx(counts_filtered, "C:\\Users\\hp\\Desktop\\R studio\\DEG_results.xlsx", rowNames = TRUE)

##########################################
### 4. PCA 
##########################################
# PCA or t-SNE Plots: To show overall sample variance and clustering based on gene expression profiles.
#Calculate PCA on the normalized data (ntd) 
pcaData <- plotPCA(ntd, intgroup = c("time", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create a PCA plot with ggplot2
ggplot(pcaData, aes(PC1, PC2, color = time, shape = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples")


#4.2 PCA or t-SNE Plots: To show overall sample variance and clustering based on gene expression profiles.
# Calculate PCA on the normalized data (ntd) --- without 72h
ntd_filtered <- ntd[, colData(dds)$time != "72h"]

pcaData_filtered <- plotPCA(ntd_filtered, intgroup = c("time", "condition"), returnData = TRUE)
percentVar_filtered <- round(100 * attr(pcaData_filtered, "percentVar"))

# Create a PCA plot with ggplot2
ggplot(pcaData_filtered, aes(PC1, PC2, color = time, shape = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Samples filtered")

#####################################
### 5.Volcano plot
#####################################
# Subset the data for 24h samples only
# Subset to 24h samples 
dds_24 <- dds[, dds$time == "24h"]

# Remove unused levels from the factors 
dds_24$time <- droplevels(dds_24$time) 
dds_24$condition <- droplevels(dds_24$condition)

# Update design to only include the 'condition' variable 
design(dds_24) <- ~ condition 

# Run DESeq2 analysis on the subset 
dds_24 <- DESeq(dds_24)

# Compare treated (LPS) vs untreated (Control) for 24h
res_24 <- results(dds_24, contrast = c("condition", "treated", "untreated"))
res_24_df <- as.data.frame(res_24)
res_24_df$gene <- rownames(res_24_df)

# Volcano plot using ggplot2
ggplot(res_24_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.4) +
  geom_point(data = subset(res_24_df, padj < 0.05 & abs(log2FoldChange) > 1),
             color = "red", alpha = 0.6) +
  ggrepel::geom_text_repel(data = subset(res_24_df, padj < 0.05 & abs(log2FoldChange) > 1),
                           aes(label = gene),
                           size = 3) +
  ggtitle("Volcano Plot for 24h Control vs LPS") +
  xlab("log2 Fold Change") +
  ylab("-log10 p-value")

######################################
# Subset the data for 48h samples only
dds_48 <- dds[, dds$time == "48h"]

# Remove unused levels from the factors 
dds_48$time <- droplevels(dds_48$time) 
dds_48$condition <- droplevels(dds_48$condition)

# Update design to only include the 'condition' variable 
design(dds_48) <- ~ condition 

# Run DESeq2 analysis on the subset 
dds_48 <- DESeq(dds_48)

# Compare treated (LPS) vs untreated (Control) for 48h
res_48 <- results(dds_48, contrast = c("condition", "treated", "untreated"))
res_48_df <- as.data.frame(res_48)
res_48_df$gene <- rownames(res_48_df)

# Volcano plot using ggplot2
ggplot(res_48_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.4) +
  geom_point(data = subset(res_48_df, padj < 0.05 & abs(log2FoldChange) > 1),
             color = "red", alpha = 0.6) +
  ggrepel::geom_text_repel(data = subset(res_48_df, padj < 0.05 & abs(log2FoldChange) > 1),
                           aes(label = gene),
                           size = 3) +
  ggtitle("Volcano Plot for 48h Control vs LPS") +
  xlab("log2 Fold Change") +
  ylab("-log10 p-value")

#########################################
# Subset the data for 72h samples only
dds_72 <- dds[, dds$time == "72h"]

# Remove unused levels from the factors 
dds_72$time <- droplevels(dds_72$time) 
dds_72$condition <- droplevels(dds_72$condition)

# Update design to only include the 'condition' variable 
design(dds_72) <- ~ condition 

# Run DESeq2 analysis on the subset 
dds_72 <- DESeq(dds_72)

# Compare treated (LPS) vs untreated (Control) for 72h
res_72 <- results(dds_72, contrast = c("condition", "treated", "untreated"))
res_72_df <- as.data.frame(res_72)
res_72_df$gene <- rownames(res_72_df)

# Volcano plot using ggplot2
ggplot(res_72_df, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(alpha = 0.4) +
  geom_point(data = subset(res_72_df, padj < 0.05 & abs(log2FoldChange) > 1),
             color = "red", alpha = 0.6) +
  ggrepel::geom_text_repel(data = subset(res_72_df, padj < 0.05 & abs(log2FoldChange) > 1),
                           aes(label = gene),
                           size = 3) +
  ggtitle("Volcano Plot for 72h Control vs LPS") +
  xlab("log2 Fold Change") +
  ylab("-log10 p-value")

#######################################
### 6. HEATMAP
#######################################
# Assume you have a vector of selected gene names
selected_genes <- c("MX1", "MX2", "ISG20", "CXCL8", "IFIT1", "IFIT3","IFI16", "NFKBIA")

# Define the desired column order
ordered_columns <- c("24h C", "24h LPS", 
                     "48h C", "48h LPS", 
                     "72h C", "72h LPS")

# Subset and reorder the normalized count data (using assay(ntd))
selected_data <- assay(ntd)[selected_genes, ordered_columns]

# Draw the heatmap
pheatmap(selected_data,
         scale = "row",          
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Heatmap of Selected Genes",
         cluster_rows = FALSE,  # Disable row clustering
         cluster_cols = FALSE)  # Disable column clustering

#######################################
### 7. GO PLOT
#######################################
# If not installed:
# BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)
# Extract significant genes from the 24h comparison
sig_genes <- rownames(subset(res_24, padj < 0.05 & abs(log2FoldChange) > 1))

# Perform GO enrichment analysis (Biological Process as an example)
ego <- enrichGO(gene          = sig_genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# Visualize the top GO terms
dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment (BP) for 24h Comparison")

##########################################
# Extract significant genes from the 48h comparison
sig_genes_48h <- rownames(subset(res_48, padj < 0.05 & abs(log2FoldChange) > 1))

# Perform GO enrichment analysis (Biological Process as an example)
ego <- enrichGO(gene          = sig_genes_48h,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# Visualize the top GO terms
dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment (BP) for 48h Comparison")

############################################
# Extract significant genes from the 72h comparison
sig_genes_72h <- rownames(subset(res_72, padj < 0.05 & abs(log2FoldChange) > 1))

# Perform GO enrichment analysis (Biological Process as an example)
ego <- enrichGO(gene          = sig_genes_72h,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.05,
                readable      = TRUE)

# Visualize the top GO terms
dotplot(ego, showCategory = 10) + ggtitle("GO Enrichment (BP) for 72h Comparison")

################################
### 8. Boxplot
################################
library(ggplot2)
library(ggpubr)
# Extract expression for one gene from the normalized data
gene_name <- "NFKBIA"
gene_expr <- assay(ntd)[gene_name, ]
gene_data <- data.frame(Expression = gene_expr,
                        Condition = samples$condition,
                        Time = samples$time)

# Create a custom order for the Time and Condition combination
gene_data$Time_Condition <- factor(interaction(gene_data$Time, gene_data$Condition), 
                                   levels = c("24h.untreated", "24h.treated", 
                                              "48h.untreated", "48h.treated", 
                                              "72h.untreated", "72h.treated"))

# Define the comparisons (untreated vs treated for each time point)
comparisons <- list(c("24h.untreated", "24h.treated"),
                    c("48h.untreated", "48h.treated"),
                    c("72h.untreated", "72h.treated"))

# Create the boxplot with significance brackets
ggplot(gene_data, aes(x = Time_Condition, y = Expression, fill = Condition)) +
  geom_boxplot() +  # Draws the boxplot
  geom_jitter(position = position_jitter(0.2), alpha = 0.4) +  # Adds scatter points for visibility
  scale_fill_manual(values = c("untreated" = "#83db81", "treated" = "#f5ef84")) +  # Custom colors
  xlab("Time and Treatment") +
  ylab("Normalized Expression") +
  ggtitle(paste("Expression of", gene_name, "across Timepoints and Treatments")) +
  theme_minimal() +  # Clean theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  
  # Add brackets with significance stars for each time point
  stat_compare_means(comparisons = comparisons, label = "p.signif", method = "t.test")  

############################################################
# Generate individual boxplots for each gene
# "untreated" vs "treated" for 24, 48 and 72 hour
comparisons <- list(c("24h.untreated", "24h.treated"),
                    c("48h.untreated", "48h.treated"),
                    c("72h.untreated", "72h.treated"))

# list of 8 genes boxplot 
boxplot_list <- lapply(selected_genes, function(gene) {
  gene_expr <- assay(ntd)[gene, ]
  gene_data <- data.frame(Expression = gene_expr, Condition = samples$condition, Time = samples$time)
  gene_data$Time_Condition <- factor(interaction(gene_data$Time, gene_data$Condition), 
                                     levels = c("24h.untreated", "24h.treated", 
                                                "48h.untreated", "48h.treated", 
                                                "72h.untreated", "72h.treated"))
  
  # Generate individual boxplots for each gene
  p <- ggplot(gene_data, aes(x = Time_Condition, y = Expression, fill = Condition)) +
    geom_boxplot() +
    geom_jitter(position = position_jitter(0.2), alpha = 0.4) +
    scale_fill_manual(values = c("untreated" = "#83db81", "treated" = "#f5ef84")) +  # Custom colors
    ggtitle(paste("Expression of", gene)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    
     # Add signs of significance with square brackets and asterisks
    stat_compare_means(comparisons = comparisons, label = "p.signif", method = "t.test", 
                       aes(group = Condition))  
  
  return(p)
})

# Arrange all boxplots in a grid
combined_plot <- gridExtra::grid.arrange(
  grobs = boxplot_list,   
  ncol = 2,              
  nrow = 4,              
  top = "Boxplot of genes of main interest"  
)

# Save image 
ggsave("Boxplot_8_genes.png", plot = combined_plot, width = 10, height = 14, dpi = 300)

##################################
### 9. venn diagram 
##################################

if (!require(VennDiagram)) install.packages("VennDiagram")
library(VennDiagram)
# Extract significant genes for each time point
sig_genes_24 <- rownames(subset(res_24, padj < 0.05 & abs(log2FoldChange) > 1))
sig_genes_48h <- rownames(subset(res_48, padj < 0.05 & abs(log2FoldChange) > 1))
sig_genes_72h <- rownames(subset(res_72, padj < 0.05 & abs(log2FoldChange) > 1))
# Create a list with the gene sets
gene_sets <- list("24h LPS" = sig_genes_24,
                  "48h LPS" = sig_genes_48h,
                  "72h LPS" = sig_genes_72h)

# Draw the Venn diagram and save it to a file
venn.diagram(
  x = gene_sets,
  filename = "C:\\Users\\hp\\Desktop\\R studio\\Venn_SignificantGenes.png", # adjust path as needed
  fill = c("#9B4F96", "#00B5B8", "#FFC107"),
  alpha = 0.50,
  cex = 1.5,
  cat.cex = 1.5,
  main = "Overlap of Significant Genes across LPS Treatments"
)


####################################
#### 10. Analysis for 24h Time Point
####################################

# Extract normalized counts from DESeq2 for the 24h time point
norm_counts_24h <- counts(dds_24, normalized = TRUE)

# Compute the mean expression for each gene in untreated (C) and treated (LPS) groups
mean_24h_C <- rowMeans(norm_counts_24h[, dds_24$condition == "untreated"], na.rm = TRUE)
mean_24h_LPS <- rowMeans(norm_counts_24h[, dds_24$condition == "treated"], na.rm = TRUE)

# Create a dataframe to store the results
df_24h <- data.frame(
  gene = rownames(norm_counts_24h),
  mean_24h_C = mean_24h_C,
  mean_24h_LPS = mean_24h_LPS,
  percentage_change = ((mean_24h_LPS - mean_24h_C) / mean_24h_C) * 100
)

# Remove infinite values caused by division by zero
df_24h <- df_24h[!is.infinite(df_24h$percentage_change), ]

# Classify genes as Increased, Decreased, or No Change based on percentage change
df_24h$Change_24h <- ifelse(df_24h$percentage_change > 0, "Increase",
                            ifelse(df_24h$percentage_change < 0, "Decrease", "No Change"))

# Keep only genes that show a significant change (Increase or Decrease)
increase_decrease_24h <- df_24h %>% filter(Change_24h %in% c("Increase", "Decrease"))

# Filter genes with a percentage change similar to STAT1 (±10%)
stat1_expression <- df_24h %>% filter(gene == "STAT1") %>% pull(percentage_change)
lower_bound <- stat1_expression * 0.9
upper_bound <- stat1_expression * 1.1
stat1_similar_genes <- increase_decrease_24h %>% 
  filter(percentage_change >= lower_bound & percentage_change <= upper_bound)

# Further filter genes based on mean expression similarity to STAT1 (±40%)
threshold <- 0.4
stat1_values <- df_24h %>% filter(gene == "STAT1") %>% dplyr::select(mean_24h_C, mean_24h_LPS)

lower_bound_C <- stat1_values$mean_24h_C * (1 - threshold)
upper_bound_C <- stat1_values$mean_24h_C * (1 + threshold)
lower_bound_LPS <- stat1_values$mean_24h_LPS * (1 - threshold)
upper_bound_LPS <- stat1_values$mean_24h_LPS * (1 + threshold)

similar_genes_24h <- increase_decrease_24h %>% 
  filter(mean_24h_C >= lower_bound_C & mean_24h_C <= upper_bound_C &
           mean_24h_LPS >= lower_bound_LPS & mean_24h_LPS <= upper_bound_LPS)

# Save the list of genes similar to STAT1
write.xlsx(similar_genes_24h, "similar_genes_24h.xlsx", rowNames = FALSE)


#### 10A. 24h Ensembl Annotation ####

# Set up the Ensembl connection and annotate the 24h gene list
similar_genes_24h_annot <- similar_genes_24h %>% 
  dplyr::rename(gene_symbol = gene)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotations_24h <- getBM(attributes = c("hgnc_symbol", "gene_biotype", "go_id", 
                                             "name_1006", "definition_1006"),
                              filters = "hgnc_symbol",
                              values = similar_genes_24h_annot$gene_symbol,
                              mart = ensembl)

# Filter annotations by pathway keywords
gene_annotations_24h <- as_tibble(gene_annotations_24h)

immune_annotations_24h <- gene_annotations_24h %>% 
  dplyr::filter(grepl("immune", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

lipid_annotations_24h <- gene_annotations_24h %>% 
  dplyr::filter(grepl("lipid", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

mitochondria_annotations_24h <- gene_annotations_24h %>% 
  dplyr::filter(grepl("mitochondria", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

combined_annotations_24h <- bind_rows(
  immune_annotations_24h,
  lipid_annotations_24h,
  mitochondria_annotations_24h
) %>% 
  dplyr::distinct()


write.xlsx(combined_annotations_24h, "combined_annotations_24h.xlsx", rowNames = FALSE)

###################################
#### 11. Analysis for 48h Time Point
###################################

# Extract normalized counts from DESeq2 for the 48h time point
norm_counts_48h <- counts(dds_48, normalized = TRUE)

# Compute the mean expression for each gene in untreated (C) and treated (LPS) groups
mean_48h_C <- rowMeans(norm_counts_48h[, dds_48$condition == "untreated"], na.rm = TRUE)
mean_48h_LPS <- rowMeans(norm_counts_48h[, dds_48$condition == "treated"], na.rm = TRUE)

# Create a dataframe to store the results
df_48h <- data.frame(
  gene = rownames(norm_counts_48h),
  mean_48h_C = mean_48h_C,
  mean_48h_LPS = mean_48h_LPS,
  percentage_change = ((mean_48h_LPS - mean_48h_C) / mean_48h_C) * 100
)

# Remove infinite values caused by division by zero
df_48h <- df_48h[!is.infinite(df_48h$percentage_change), ]

# Classify genes as Increased, Decreased, or No Change based on percentage change
df_48h$Change_48h <- ifelse(df_48h$percentage_change > 0, "Increase",
                            ifelse(df_48h$percentage_change < 0, "Decrease", "No Change"))

# Keep only genes that show a significant change (Increase or Decrease)
increase_decrease_48h <- df_48h %>% filter(Change_48h %in% c("Increase", "Decrease"))

# Filter genes with a percentage change similar to STAT1 (±10%)
stat1_expression <- df_48h %>% filter(gene == "STAT1") %>% pull(percentage_change)
lower_bound <- stat1_expression * 0.9
upper_bound <- stat1_expression * 1.1
stat1_similar_genes <- increase_decrease_48h %>% 
  filter(percentage_change >= lower_bound & percentage_change <= upper_bound)

# Further filter genes based on mean expression similarity to STAT1 (±40%)
threshold <- 0.4
stat1_values <- df_48h %>% filter(gene == "STAT1") %>% dplyr::select(mean_48h_C, mean_48h_LPS)

lower_bound_C <- stat1_values$mean_48h_C * (1 - threshold)
upper_bound_C <- stat1_values$mean_48h_C * (1 + threshold)
lower_bound_LPS <- stat1_values$mean_48h_LPS * (1 - threshold)
upper_bound_LPS <- stat1_values$mean_48h_LPS * (1 + threshold)

similar_genes_48h <- increase_decrease_48h %>% 
  filter(mean_48h_C >= lower_bound_C & mean_48h_C <= upper_bound_C &
           mean_48h_LPS >= lower_bound_LPS & mean_48h_LPS <= upper_bound_LPS)

# Save the list of genes similar to STAT1
write.xlsx(similar_genes_48h, "similar_genes_48h.xlsx", rowNames = FALSE)

#### 11A. 48h Ensembl Annotation ####

# Set up the Ensembl connection and annotate the 24h gene list
similar_genes_48h_annot <- similar_genes_48h %>% 
  dplyr::rename(gene_symbol = gene)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotations_48h <- getBM(attributes = c("hgnc_symbol", "gene_biotype", "go_id", 
                                             "name_1006", "definition_1006"),
                              filters = "hgnc_symbol",
                              values = similar_genes_48h_annot$gene_symbol,
                              mart = ensembl)

# Filter annotations by pathway keywords
gene_annotations_48h <- as_tibble(gene_annotations_48h)

immune_annotations_48h <- gene_annotations_48h %>% 
  dplyr::filter(grepl("immune", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

lipid_annotations_48h <- gene_annotations_48h %>% 
  dplyr::filter(grepl("lipid", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

mitochondria_annotations_48h <- gene_annotations_48h %>% 
  dplyr::filter(grepl("mitochondria", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

combined_annotations_48h <- bind_rows(
  immune_annotations_48h,
  lipid_annotations_48h,
  mitochondria_annotations_48h
) %>% 
  dplyr::distinct()


write.xlsx(combined_annotations_48h, "combined_annotations_48h.xlsx", rowNames = FALSE)

###################################
#### 12. Analysis for 72h Time Point
###################################

# Extract normalized counts from DESeq2 for the 48h time point
norm_counts_72h <- counts(dds_72, normalized = TRUE)

# Compute the mean expression for each gene in untreated (C) and treated (LPS) groups
mean_72h_C <- rowMeans(norm_counts_72h[, dds_72$condition == "untreated"], na.rm = TRUE)
mean_72h_LPS <- rowMeans(norm_counts_72h[, dds_72$condition == "treated"], na.rm = TRUE)

# Create a dataframe to store the results
df_72h <- data.frame(
  gene = rownames(norm_counts_72h),
  mean_72h_C = mean_72h_C,
  mean_72h_LPS = mean_72h_LPS,
  percentage_change = ((mean_72h_LPS - mean_72h_C) / mean_72h_C) * 100
)

# Remove infinite values caused by division by zero
df_72h <- df_72h[!is.infinite(df_72h$percentage_change), ]

# Classify genes as Increased, Decreased, or No Change based on percentage change
df_72h$Change_72h <- ifelse(df_72h$percentage_change > 0, "Increase",
                            ifelse(df_72h$percentage_change < 0, "Decrease", "No Change"))

# Keep only genes that show a significant change (Increase or Decrease)
increase_decrease_72h <- df_72h %>% filter(Change_72h %in% c("Increase", "Decrease"))

# Filter genes with a percentage change similar to STAT1 (±10%)
stat1_expression <- df_72h %>% filter(gene == "STAT1") %>% pull(percentage_change)
lower_bound <- stat1_expression * 0.9
upper_bound <- stat1_expression * 1.1
stat1_similar_genes <- increase_decrease_72h %>% 
  filter(percentage_change >= lower_bound & percentage_change <= upper_bound)

# Further filter genes based on mean expression similarity to STAT1 (±40%)
threshold <- 0.4
stat1_values <- df_72h %>% filter(gene == "STAT1") %>% dplyr::select(mean_72h_C, mean_72h_LPS)

lower_bound_C <- stat1_values$mean_72h_C * (1 - threshold)
upper_bound_C <- stat1_values$mean_72h_C * (1 + threshold)
lower_bound_LPS <- stat1_values$mean_72h_LPS * (1 - threshold)
upper_bound_LPS <- stat1_values$mean_72h_LPS * (1 + threshold)

similar_genes_72h <- increase_decrease_72h %>% 
  filter(mean_72h_C >= lower_bound_C & mean_72h_C <= upper_bound_C &
           mean_72h_LPS >= lower_bound_LPS & mean_72h_LPS <= upper_bound_LPS)

# Save the list of genes similar to STAT1
write.xlsx(similar_genes_72h, "similar_genes_72h.xlsx", rowNames = FALSE)

#### 12A. 72h Ensembl Annotation ####

# Set up the Ensembl connection and annotate the 72h gene list
similar_genes_72h_annot <- similar_genes_72h %>% 
  dplyr::rename(gene_symbol = gene)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

gene_annotations_72h <- getBM(attributes = c("hgnc_symbol", "gene_biotype", "go_id", 
                                             "name_1006", "definition_1006"),
                              filters = "hgnc_symbol",
                              values = similar_genes_72h_annot$gene_symbol,
                              mart = ensembl)

# Filter annotations by pathway keywords
gene_annotations_72h <- as_tibble(gene_annotations_72h)

immune_annotations_72h <- gene_annotations_72h %>% 
  dplyr::filter(grepl("immune", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

lipid_annotations_72h <- gene_annotations_72h %>% 
  dplyr::filter(grepl("lipid", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

mitochondria_annotations_72h <- gene_annotations_72h %>% 
  dplyr::filter(grepl("mitochondria", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>% 
  dplyr::distinct()

combined_annotations_72h <- bind_rows(
  immune_annotations_72h,
  lipid_annotations_72h,
  mitochondria_annotations_72h
) %>% 
  dplyr::distinct()


write.xlsx(combined_annotations_72h, "combined_annotations_72h.xlsx", rowNames = FALSE)


##############################################################################
####### top 10 upregulated/downregulated genes at each time point
##############################################################################
library(ggpubr)
####################### 24h ##################################################
# ---- Top 10 UP ----
top_up_24 <- res_24_df %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10) %>%
  select(Gene = gene, `Log2FC` = log2FoldChange, `Adj.Pval` = padj)

tab_up_24 <- ggtexttable(top_up_24,
                         theme = ttheme("light"),
                         rows = NULL) %>%
  tab_add_title(text = "Top 10 Upregulated Genes (24h)", face = "bold")

# ---- Top 10 DOWN ----
top_down_24 <- res_24_df %>%
  filter(padj < 0.05, log2FoldChange < -1) %>%
  arrange(log2FoldChange) %>%
  head(10) %>%
  select(Gene = gene, `Log2FC` = log2FoldChange, `Adj.Pval` = padj)

if (nrow(top_down_24) > 0) {
  tab_down_24 <- ggtexttable(top_down_24,
                             theme = ttheme("light"),
                             rows = NULL) %>%
    tab_add_title(text = "Top 10 Downregulated Genes (24h)", face = "bold")
} else {
  tab_down_24 <- ggtexttable(data.frame(Note = "No significant downregulated genes"),
                             theme = ttheme("light")) %>%
    tab_add_title(text = "Top 10 Downregulated Genes (24h)", face = "bold")
}

# ---- Combine and Save ----
library(grid)
png("Top_24h_DEGs_Tables.png", width = 1200, height = 600, res = 150)
grid.draw(plot_24_tables)
dev.off()

####################### 48h ##################################################
# ---- Top 10 UP ----
top_up_48 <- res_48_df %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10) %>%
  select(Gene = gene, `Log2FC` = log2FoldChange, `Adj.Pval` = padj)

tab_up_48 <- ggtexttable(top_up_48,
                         theme = ttheme("light"),
                         rows = NULL) %>%
  tab_add_title(text = "Top 10 Upregulated Genes (48h)", face = "bold")

# ---- Top 10 DOWN ----
top_down_48 <- res_48_df %>%
  filter(padj < 0.05, log2FoldChange < -1) %>%
  arrange(log2FoldChange) %>%
  head(10) %>%
  select(Gene = gene, `Log2FC` = log2FoldChange, `Adj.Pval` = padj)

if (nrow(top_down_48) > 0) {
  tab_down_48 <- ggtexttable(top_down_48,
                             theme = ttheme("light"),
                             rows = NULL) %>%
    tab_add_title(text = "Top 10 Downregulated Genes (48h)", face = "bold")
} else {
  tab_down_48 <- ggtexttable(data.frame(Note = "No significant downregulated genes"),
                             theme = ttheme("light")) %>%
    tab_add_title(text = "Top 10 Downregulated Genes (48h)", face = "bold")
}

# ---- Combine and Save ----
png("Top_48h_DEGs_Tables.png", width = 1200, height = 600, res = 150)
grid.draw(plot_48_tables)
dev.off()

####################### 72h ##################################################
# ---- Top 10 UP ----
top_up_72 <- res_72_df %>%
  filter(padj < 0.05, log2FoldChange > 1) %>%
  arrange(desc(log2FoldChange)) %>%
  head(10) %>%
  select(Gene = gene, `Log2FC` = log2FoldChange, `Adj.Pval` = padj)

tab_up_72 <- ggtexttable(top_up_72,
                         theme = ttheme("light"),
                         rows = NULL) %>%
  tab_add_title(text = "Top 10 Upregulated Genes (72h)", face = "bold")

# ---- Top 10 DOWN ----
top_down_72 <- res_72_df %>%
  filter(padj < 0.05, log2FoldChange < -1) %>%
  arrange(log2FoldChange) %>%
  head(10) %>%
  select(Gene = gene, `Log2FC` = log2FoldChange, `Adj.Pval` = padj)

if (nrow(top_down_72) > 0) {
  tab_down_72 <- ggtexttable(top_down_72,
                             theme = ttheme("light"),
                             rows = NULL) %>%
    tab_add_title(text = "Top 10 Downregulated Genes (72h)", face = "bold")
} else {
  tab_down_72 <- ggtexttable(data.frame(Note = "No significant downregulated genes"),
                             theme = ttheme("light")) %>%
    tab_add_title(text = "Top 10 Downregulated Genes (72h)", face = "bold")
}

# ---- Combine and Save ----
png("Top_72h_DEGs_Tables.png", width = 1200, height = 600, res = 150)
grid.draw(plot_72_tables)
dev.off()




