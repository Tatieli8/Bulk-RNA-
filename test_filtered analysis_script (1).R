#### to get to 'Thomas addition'####

####Loading Libraries####

#Requirements#
if (!require("BiocManager", quietly = TRUE)) + install.packages("BiocManager")
BiocManager::install(version = "3.18")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("pheatmap")
BiocManager::install("vsn")
install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("readr")
install.packages("biomaRt")

#library loading
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)
library(DESeq2)
library("ggplot2")
library("pheatmap")
library("vsn")
library(RColorBrewer)
library(readr)

####Loading files####

Kimdata <- read.csv("C:\\Users\\hp\\Desktop\\DANDRITE\\Script\\Thomas.rawcounts.csv")

#S54 = 48h LPS
#S114 = 48h C
#S125 = 72h LPS
#S12 = 48h C Rep1
#S36 = 24h C
#S122 = 72h C
#S24 = 24h C Rep1
#S108 = 24h LPS
#S123 = 72h C Rep1
#S18 = 48h LPS Rep1
#S107 = 24h C Rep2
#S42 = 24h LPS Rep1
#S30 = 24h LPS Rep2
#S124 = 72h LPS Rep1
#S116 = 48h LPS Rep2
#S48 = 48h C Rep2

treatment <- c("X", "48h LPS", "48h C", "72h LPS", "48h C Rep1",
               "24h C", "72h C", "24h C Rep1",
               "24h LPS", "72h C Rep1", "48h LPS Rep1",
               "24h C Rep2", "24h LPS Rep1", "24h LPS Rep2",
               "72h LPS Rep1", "48h LPS Rep2", "48h C Rep2" )
colnames(Kimdata) <- treatment
# Reorder columns
Kimdata <- Kimdata[, c("X", "24h C", "24h C Rep1", "24h C Rep2", 
                       "24h LPS", "24h LPS Rep1", "24h LPS Rep2",
                       "48h C", "48h C Rep1", "48h C Rep2", 
                       "48h LPS", "48h LPS Rep1", "48h LPS Rep2",
                       "72h C", "72h C Rep1", "72h LPS", "72h LPS Rep1")]
#X
#24h C = untreated
#24h C Rep1 = untreated
#24h C Rep2 = untreated
#24h LPS = treated
#24h LPS Rep1 = treated
#24h LPS Rep2 = treated
#48h C = untreated
#48h C Rep1 = untreated
#48h C Rep2 = untreated
#48h LPS = treated
#48h LPS Rep1 = treated
#48h LPS Rep2 = treated
#72h C = untreated
#72h C Rep1 = untreated
#72h LPS = treated
#72h LPS Rep1 = treated

###Converting Ensembl Codes to gene names using org.Hs.eg.db###
library(dplyr)
library(org.Hs.eg.db)
library(biomaRt)


geneID <- select(org.Hs.eg.db, keys = Kimdata$X, keytype = "ENSEMBL", columns = "SYMBOL")
geneID2 <- geneID[!duplicated(geneID[c('ENSEMBL')]), ]
geneID3 <- na.omit(geneID2)
Symbols_old <- geneID2$SYMBOL # contains NA values
Symbols <- geneID3$SYMBOL # without NA values
test <- subset(Kimdata, select = -c(X))
Rows <- dplyr::coalesce(Symbols_old, rownames(test))    
Unique.Rows <- make.unique(Rows, sep = " ")
rownames(test) <- Unique.Rows
test <- test[rownames(test) %in% Symbols,]

#Conditions & Treatment Time
condition <- c("untreated", "untreated", "untreated",
               "treated", "treated", "treated",
               "untreated", "untreated", "untreated",
               "treated", "treated", "treated", 
               "untreated", "untreated", "treated", "treated")
time <- c("24h", "24h", "24h","24h", "24h", "24h", 
          "48h", "48h", "48h", "48h", "48h", "48h", 
          "72h", "72h","72h", "72h")
samples <- cbind(condition, time)
rownames(samples) <- colnames(test)
samples <- data.frame(samples)

### Filtering & DESeq2 
library(DESeq2)
test_filtered <- test %>%
  filter(rowSums(across(where(is.numeric))) != 0) 

dds <- DESeqDataSetFromMatrix(countData = test_filtered,
                              colData = samples,
                              design= ~ condition + time)
dds <- DESeq(dds)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

ntd <- normTransform(dds)
df <- as.data.frame(colData(dds)[,c("condition","time")])


####export in a excel file#####

install.packages("xlsx")
install.packages("openxlsx")
library(xlsx)
library(openxlsx)

install.packages("writexl")
library(writexl)
write.xlsx(test_filtered, "test_filtered.xlsx", rowNames=TRUE)

####changes in tripletes in 24h####
# Load the necessary library
library(dplyr)

# Check the column names in the dataset
names(test_filtered)

# Calculate the mean for 24h C and 24h LPS
test_filtered$mean_24h_C <- rowMeans(test_filtered[, 
                       c('24h C', '24h C Rep1', '24h C Rep2')], na.rm = TRUE)
test_filtered$mean_24h_LPS <- rowMeans(test_filtered[, 
                      c('24h LPS', '24h LPS Rep1', '24h LPS Rep2')], na.rm = TRUE)

# Calculate the percentage change from 24h C to 24h LPS
test_filtered$percentage_change <- ((test_filtered$mean_24h_LPS - test_filtered$mean_24h_C) / test_filtered$mean_24h_C) * 100

# Identify if the change is an increase or a decrease
test_filtered$Change_24h <- ifelse(test_filtered$percentage_change > 0, "Increase",
                                   ifelse(test_filtered$percentage_change < 0, "Decrease", "No Change"))

# Remove rows with Inf or -Inf in percentage_change
test_filtered <- test_filtered[!is.infinite(test_filtered$percentage_change), ]

# Keep only the rows where the gene increased or decreased
increase_decrease_24h <- test_filtered[test_filtered$Change_24h %in% c("Increase", "Decrease"), ]
increase_decrease_24h <- increase_decrease_24h[, -c(7, 8, 9, 10, 11, 12, 13, 14, 15, 16)]  # delete others columns 

# Check the structure of the dataset to confirm it's a data frame
head(increase_decrease_24h)
str(increase_decrease_24h)
colnames(increase_decrease_24h)


# Get STAT1's average expression 
stat1_expression <- increase_decrease_24h %>%
  dplyr::filter(rownames(increase_decrease_24h) == "STAT1")  %>%  # Use dplyr::filter to avoid conflicts
  pull(percentage_change)

# Set the proximity threshold (for example, within 10% of STAT1's expression)
lower_bound <- stat1_expression * 0.9
upper_bound <- stat1_expression * 1.1

# Filter genes with average expression near STAT1
stat1_similargene <- increase_decrease_24h %>%
  dplyr::filter(percentage_change >= lower_bound & percentage_change <= upper_bound)


# Get STAT1's mean values
# Load necessary libraries
library(readr)
library(dplyr)

stat1_values <- increase_decrease_24h %>%
  dplyr::filter(rownames(increase_decrease_24h) == "STAT1")  %>%
  dplyr::select(mean_24h_C, mean_24h_LPS)  # Extract both mean values

# Store the values for easy access
mean_24h_C_stat1 <- stat1_values$mean_24h_C
mean_24h_LPS_stat1 <- stat1_values$mean_24h_LPS

# Set proximity thresholds (for example, within 10% of STAT1's values)
threshold_C <- 0.1  # 10% threshold
threshold_LPS <- 0.1  # 10% threshold

lower_bound_C <- mean_24h_C_stat1 * (1 - threshold_C)
upper_bound_C <- mean_24h_C_stat1 * (1 + threshold_C)

lower_bound_LPS <- mean_24h_LPS_stat1 * (1 - threshold_LPS)
upper_bound_LPS <- mean_24h_LPS_stat1 * (1 + threshold_LPS)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_24h <- increase_decrease_24h %>%
  dplyr::filter(mean_24h_C >= lower_bound_C & mean_24h_C <= upper_bound_C &
                  mean_24h_LPS >= lower_bound_LPS & mean_24h_LPS <= upper_bound_LPS)


#need to amplify the range 
# Set proximity thresholds (for example, within 30% of STAT1's values)
threshold_C <- 0.3  # 30% threshold
threshold_LPS <- 0.3  # 30% threshold

lower_bound_C <- mean_24h_C_stat1 * (1 - threshold_C)
upper_bound_C <- mean_24h_C_stat1 * (1 + threshold_C)

lower_bound_LPS <- mean_24h_LPS_stat1 * (1 - threshold_LPS)
upper_bound_LPS <- mean_24h_LPS_stat1 * (1 + threshold_LPS)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_24h <- increase_decrease_24h %>%
  dplyr::filter(mean_24h_C >= lower_bound_C & mean_24h_C <= upper_bound_C &
                  mean_24h_LPS >= lower_bound_LPS & mean_24h_LPS <= upper_bound_LPS)


write.xlsx(similar_genes, "similar_genes_excel.xlsx", rowNames=TRUE)


#need to amplify the range 
# Set proximity thresholds (for example, within 40% of STAT1's values)
threshold_C <- 0.4  # 40% threshold
threshold_LPS <- 0.4  # 40% threshold

lower_bound_C <- mean_24h_C_stat1 * (1 - threshold_C)
upper_bound_C <- mean_24h_C_stat1 * (1 + threshold_C)

lower_bound_LPS <- mean_24h_LPS_stat1 * (1 - threshold_LPS)
upper_bound_LPS <- mean_24h_LPS_stat1 * (1 + threshold_LPS)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_24h <- increase_decrease_24h %>%
  dplyr::filter(mean_24h_C >= lower_bound_C & mean_24h_C <= upper_bound_C &
                  mean_24h_LPS >= lower_bound_LPS & mean_24h_LPS <= upper_bound_LPS)


write.xlsx(similar_genes_24h, "similar_genes_24h.xlsx", rowNames=TRUE)

####24h Ensembl####

library(biomaRt)
library(dplyr)
library(openxlsx)
library(tibble)

#set up Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

similar_genes_24h <- rownames_to_column(similar_genes_24h, var = "gene_symbol")

#extract pathway and annotation information for genes of interest (24h)
gene_annotations_24h <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype", "go_id", "name_1006", "definition_1006"),
  filters = "hgnc_symbol",
  values = similar_genes_24h, #replace the data with your data set 
  mart = ensembl
)

# Filter annotations for specific pathways
immune_annotations_24h <- gene_annotations_24h %>% 
  filter(grepl("immune", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>%
  distinct()

lipid_annotations_24h <- gene_annotations_24h %>% 
  filter(grepl("lipid", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>%
  distinct()

mitochondria_annotations_24h <- gene_annotations_24h %>% 
  filter(grepl("mitochondria", name_1006, ignore.case = TRUE)) %>%
  dplyr::select(hgnc_symbol, gene_biotype, name_1006, definition_1006) %>%
  distinct()

# Combine all annotations into a single data set
combined_annotations_24h <- bind_rows(
  immune_annotations_24h,
  lipid_annotations_24h,
  mitochondria_annotations_24h
) %>%
  distinct()  # Remove duplicate rows

# Save to Excel
write.xlsx(combined_annotations_24h, "combined_annotations_24h.xlsx", rowNames = FALSE)


####changes in tripletes in 48h####
# Load the necessary library
library(dplyr)

# Calculate the mean for 48h C and 48h LPS
test_filtered$mean_48h_C <- rowMeans(test_filtered[, 
                                                   c('48h C', '48h C Rep1', '48h C Rep2')], na.rm = TRUE)
test_filtered$mean_48h_LPS <- rowMeans(test_filtered[, 
                                                     c('48h LPS', '48h LPS Rep1', '48h LPS Rep2')], na.rm = TRUE)

# Calculate the percentage change from 48h C to 48h LPS
test_filtered$percentage_change <- ((test_filtered$mean_48h_LPS - test_filtered$mean_48h_C) / test_filtered$mean_48h_C) * 100

# Identify if the change is an increase or a decrease
test_filtered$Change_48h <- ifelse(test_filtered$percentage_change > 0, "Increase",
                                   ifelse(test_filtered$percentage_change < 0, "Decrease", "No Change"))

# Remove rows with Inf or -Inf in percentage_change
test_filtered <- test_filtered[!is.infinite(test_filtered$percentage_change), ]

# Keep only the rows where the gene increased or decreased
increase_decrease_48h <- test_filtered[test_filtered$Change_48h %in% c("Increase", "Decrease"), ]
increase_decrease_48h <- increase_decrease_48h[, -c(1, 2, 3, 4, 5, 6, 13, 14, 15, 16)]  # delete others columns 

# Check the structure of the dataset to confirm it's a data frame
str(increase_decrease_48h)
colnames(increase_decrease_48h)
head(increase_decrease_48h)

# Get STAT1's average expression 
stat1_expression_48 <- increase_decrease_48h %>%
  dplyr::filter(rownames(increase_decrease_48h) == "STAT1") %>%  # Use dplyr::filter to avoid conflicts
  pull(percentage_change)

# Set the proximity threshold (for example, within 10% of STAT1's expression)
lower_bound_48 <- stat1_expression_48 * 0.9
upper_bound_48 <- stat1_expression_48 * 1.1

# Filter genes with average expression near STAT1
similar_genes_48 <- increase_decrease_48h %>%
  dplyr::filter(percentage_change >= lower_bound_48 & percentage_change <= upper_bound_48)


# Get STAT1's mean values at 48h
# Load necessary libraries
library(readr)
library(dplyr)

stat1_values_48 <- increase_decrease_48h %>%
  dplyr::filter(rownames(increase_decrease_48h) == "STAT1") %>%
  dplyr::select(mean_48h_C, mean_48h_LPS)  # Extract both mean values

# Store the values for easy access
mean_48h_C_stat1 <- stat1_values_48$mean_48h_C
mean_48h_LPS_stat1 <- stat1_values_48$mean_48h_LPS

# Set proximity thresholds (for example, within 10% of STAT1's values)
threshold_C_48 <- 0.1  # 10% threshold
threshold_LPS_48 <- 0.1  # 10% threshold

lower_bound_C_48 <- mean_48h_C_stat1 * (1 - threshold_C_48)
upper_bound_C_48 <- mean_48h_C_stat1 * (1 + threshold_C_48)

lower_bound_LPS_48 <- mean_48h_LPS_stat1 * (1 - threshold_LPS_48)
upper_bound_LPS_48 <- mean_48h_LPS_stat1 * (1 + threshold_LPS_48)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_48 <- increase_decrease_48h %>%
  dplyr::filter(mean_48h_C >= lower_bound_C_48 & mean_48h_C <= upper_bound_C_48 &
                  mean_48h_LPS >= lower_bound_LPS_48 & mean_48h_LPS <= upper_bound_LPS_48)


#need to amplify the range 
# Set proximity thresholds (for example, within 30% of STAT1's values)
threshold_C_48_0.3 <- 0.3  # 30% threshold
threshold_LPS_48_0.3 <- 0.3  # 30% threshold

lower_bound_C_48_0.3 <- mean_48h_C_stat1 * (1 - threshold_C_48_0.3)
upper_bound_C_48_0.3 <- mean_48h_C_stat1 * (1 + threshold_C_48_0.3)

lower_bound_LPS_48_0.3 <- mean_48h_LPS_stat1 * (1 - threshold_LPS_48_0.3)
upper_bound_LPS_48_0.3 <- mean_48h_LPS_stat1 * (1 + threshold_LPS_48_0.3)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_48 <- increase_decrease_48h %>%
  dplyr::filter(mean_48h_C >= lower_bound_C_48_0.3 & mean_48h_C <= upper_bound_C_48_0.3 &
                  mean_48h_LPS >= lower_bound_LPS_48_0.3 & mean_48h_LPS <= upper_bound_LPS_48_0.3)


write.xlsx(similar_genes48, "similar_genes48_excel.xlsx", rowNames=TRUE)

#need to amplify the range 40%
# Set proximity thresholds (for example, within 40% of STAT1's values)
threshold_C_48_0.4 <- 0.4  # 40% threshold
threshold_LPS_48_0.4 <- 0.4  # 40% threshold

lower_bound_C_48_0.4 <- mean_48h_C_stat1 * (1 - threshold_C_48_0.4)
upper_bound_C_48_0.4 <- mean_48h_C_stat1 * (1 + threshold_C_48_0.4)

lower_bound_LPS_48_0.4 <- mean_48h_LPS_stat1 * (1 - threshold_LPS_48_0.4)
upper_bound_LPS_48_0.4 <- mean_48h_LPS_stat1 * (1 + threshold_LPS_48_0.4)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_48h <- increase_decrease_48h %>%
  dplyr::filter(mean_48h_C >= lower_bound_C_48_0.4 & mean_48h_C <= upper_bound_C_48_0.4 &
                  mean_48h_LPS >= lower_bound_LPS_48_0.4 & mean_48h_LPS <= upper_bound_LPS_48_0.4)


write.xlsx(similar_genes_48h, "similar_genes_48hl.xlsx", rowNames=TRUE)


####48h Ensembl####
library(biomaRt)
library(dplyr)
library(openxlsx)
library(tibble)

#set up Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
similar_genes_48h <- rownames_to_column(similar_genes_48h, var = "gene_symbol")

#extract pathway and annotation information for genes of interest (48h)
gene_annotations_48h <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype", "go_id", "name_1006", "definition_1006"),
  filters = "hgnc_symbol",
  values = similar_genes_48h, #sostituisci "gene_column" con il nome della colonna dei geni
  mart = ensembl
)

# Filter annotations for specific pathways
immune_annotations_48h <- gene_annotations_48h[grepl
              ("immune", gene_annotations_48h$name_1006, ignore.case = TRUE), 
             c("hgnc_symbol", "gene_biotype", "name_1006", "definition_1006")]
immune_annotations_48h <- unique(immune_annotations_48h)


lipid_annotations_48h <- gene_annotations_48h[grepl
              ("lipid", gene_annotations_48h$name_1006, ignore.case = TRUE), 
             c("hgnc_symbol", "gene_biotype", "name_1006", "definition_1006")]
lipid_annotations_48h <- unique(lipid_annotations_48h)
 

mitochondria_annotations_48h <- gene_annotations_48h[grepl 
              ("mitochondria", gene_annotations_48h$name_1006, ignore.case = TRUE),
             c("hgnc_symbol", "gene_biotype", "name_1006", "definition_1006")]
mitochondria_annotations_48h <- unique(mitochondria_annotations_48h)

# Combine all annotations into a single data set
combined_annotations_48h <- bind_rows(
  immune_annotations_48h,
  lipid_annotations_48h,
  mitochondria_annotations_48h
) %>%
  distinct()  # Remove duplicate rows

# Save to Excel
write.xlsx(combined_annotations_48h, "ensembl_48h.xlsx", rowNames = FALSE)

install.packages("rJava")
install.packages("openxlsx")
library(xlsx)
library(openxlsx)
remove.packages("xlsx")
remove.packages("rJava")
install.packages("openxlsx")
library(openxlsx)

####create pheatmap####
gene_names <- c("MRPS14", "IFI16", "NDUFC1", "DBNL", "SDCBP",
                "PLD3", "TOMM40", "RCN3", "ZNFX1", "VAT1",
                "SEC14L1", "HEXA", "NFKBIA", "SLC43A3", "PMPCA")
# Ensure these gene names are in the row names of 'dds'
if(!all(gene_names %in% rownames(dds))) {stop("Some specified genes are not in the dataset")}


# Normalization step: Normalize the rows by Z-score
normalized_data <- t(scale(t(assay(ntd)[gene_names, ])))


# Reorder columns based on time and condition
#re-visit "sample"
condition <- c("treated", "untreated", "untreated",
               "untreated", "untreated", "treated",
               "treated", "untreated", "treated",
               "treated", "treated", "untreated")
time <- c("48h", "48h","48h", "24h", "24h", "24h",
          "48h", "24h", "24h", "24h", "48h", "48h")
sample_names <- c("48h LPS", "48h C", "48h C Rep1",
                  "24h C", "24h C Rep1",
                  "24h LPS", "48h LPS Rep1",
                  "24h C Rep2", "24h LPS Rep1", "24h LPS Rep2",
                  "48h LPS Rep2", "48h C Rep2")
samples <- data.frame(condition, time, row.names = sample_names)


#order stuff
custom_order <- c("24h C", "24h C Rep1", "24h C Rep2",
                  "24h LPS", "24h LPS Rep1", "24h LPS Rep2",
                  "48h C", "48h C Rep1", "48h C Rep2",
                  "48h LPS", "48h LPS Rep1", "48h LPS Rep2")
samples_ordered <- samples[custom_order, ]


#order data
normalized_data_ordered <- normalized_data[, match(custom_order, colnames(normalized_data))]


# Define the colors for the annotations
ann_colors <- list(
  time = c("24h" = "palegreen",
           "48h" = "gold"),
  condition = c("treated" = "cornflowerblue",
                "untreated" = "plum"))


# Create a color palette transitioning from blue to orange to red
color_palette <- colorRampPalette(rev(brewer.pal(10, "RdYlBu")))(256)
install.packages("brewer.pal")
install.packages("pheatmap")
library(pheatmap)


# Create the heatmap with normalized data and the custom ordered columns
pheatmap(normalized_data_ordered, border_color = "white",
         color = color_palette,
         cluster_rows = FALSE, show_rownames = TRUE,
         labels_row = c("MRPS14", "IFI16", "NDUFC1", "DBNL", "SDCBP",
                        "PLD3", "TOMM40", "RCN3", "ZNFX1", "VAT1",
                        "SEC14L1", "HEXA", "NFKBIA", "SLC43A3", "PMPCA"),
         cluster_cols = FALSE, annotation_col = samples_ordered,
         annotation_colors = ann_colors, angle_col = 45)

####changes in tripletes in 72h####
# Load the necessary library
library(dplyr)

# Calculate the mean for 48h C and 48h LPS
test_filtered$mean_72h_C <- rowMeans(test_filtered[, 
                                                   c('72h C', '72h C Rep1')], na.rm = TRUE)
test_filtered$mean_72h_LPS <- rowMeans(test_filtered[, 
                                                     c('72h LPS', '72h LPS Rep1')], na.rm = TRUE)

# Calculate the percentage change from 48h C to 48h LPS
test_filtered$percentage_change <- ((test_filtered$mean_72h_LPS - test_filtered$mean_72h_C) / test_filtered$mean_72h_C) * 100

# Identify if the change is an increase or a decrease
test_filtered$Change_72h <- ifelse(test_filtered$percentage_change > 0, "Increase",
                                   ifelse(test_filtered$percentage_change < 0, "Decrease", "No Change"))

# Remove rows with Inf or -Inf in percentage_change
test_filtered <- test_filtered[!is.infinite(test_filtered$percentage_change), ]

# Keep only the rows where the gene increased or decreased
increase_decrease_72h <- test_filtered[test_filtered$Change_72h %in% c("Increase", "Decrease"), ]
increase_decrease_72h <- increase_decrease_72h[, -c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)]  # delete others columns 

# Check the structure of the dataset to confirm it's a data frame
str(increase_decrease_72h)
colnames(increase_decrease_72h)
head(increase_decrease_72h)

# Get STAT1's average expression 
stat1_expression_72 <- increase_decrease_72h %>%
  dplyr::filter(rownames(increase_decrease_72h) == "STAT1") %>%  # Use dplyr::filter to avoid conflicts
  pull(percentage_change)

# Set the proximity threshold (for example, within 10% of STAT1's expression)
lower_bound_72 <- stat1_expression_72 * 0.9
upper_bound_72 <- stat1_expression_72 * 1.1

# Filter genes with average expression near STAT1
similar_genes_72 <- increase_decrease_72h %>%
  dplyr::filter(percentage_change >= lower_bound_72 & percentage_change <= upper_bound_72)


# Get STAT1's mean values at 72h
# Load necessary libraries
library(readr)
library(dplyr)

stat1_values_72 <- increase_decrease_72h %>%
  dplyr::filter(rownames(increase_decrease_72h) == "STAT1") %>%
  dplyr::select(mean_72h_C, mean_72h_LPS)  # Extract both mean values

# Store the values for easy access
mean_72h_C_stat1 <- stat1_values_72$mean_72h_C
mean_72h_LPS_stat1 <- stat1_values_72$mean_72h_LPS

# Set proximity thresholds (for example, within 10% of STAT1's values)
threshold_C_72 <- 0.1  # 10% threshold
threshold_LPS_72 <- 0.1  # 10% threshold

lower_bound_C_72 <- mean_72h_C_stat1 * (1 - threshold_C_72)
upper_bound_C_72 <- mean_72h_C_stat1 * (1 + threshold_C_72)

lower_bound_LPS_72 <- mean_72h_LPS_stat1 * (1 - threshold_LPS_72)
upper_bound_LPS_72 <- mean_72h_LPS_stat1 * (1 + threshold_LPS_72)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_72 <- increase_decrease_72h %>%
  dplyr::filter(mean_72h_C >= lower_bound_C_72 & mean_72h_C <= upper_bound_C_72 &
                  mean_72h_LPS >= lower_bound_LPS_72 & mean_72h_LPS <= upper_bound_LPS_72)


#need to amplify the range 
# Set proximity thresholds (for example, within 30% of STAT1's values)
threshold_C_72_0.3 <- 0.3  # 30% threshold
threshold_LPS_72_0.3 <- 0.3  # 30% threshold

lower_bound_C_72_0.3 <- mean_72h_C_stat1 * (1 - threshold_C_72_0.3)
upper_bound_C_72_0.3 <- mean_72h_C_stat1 * (1 + threshold_C_72_0.3)

lower_bound_LPS_72_0.3 <- mean_72h_LPS_stat1 * (1 - threshold_LPS_72_0.3)
upper_bound_LPS_72_0.3 <- mean_72h_LPS_stat1 * (1 + threshold_LPS_72_0.3)

# Filter genes with mean_24h_C and mean_24h_LPS near STAT1
similar_genes_72h <- increase_decrease_72h %>%
  dplyr::filter(mean_72h_C >= lower_bound_C_72_0.3 & mean_72h_C <= upper_bound_C_72_0.3 &
                  mean_72h_LPS >= lower_bound_LPS_72_0.3 & mean_72h_LPS <= upper_bound_LPS_72_0.3)


write.xlsx(similar_genes_72h, "similar_genes_72h.xlsx", rowNames=TRUE)

####72h Ensembl####
library(biomaRt)
library(dplyr)
library(openxlsx)
library(tibble)

#set up Ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
similar_genes_72h <- rownames_to_column(similar_genes_72h, var = "gene_symbol")

#extract pathway and annotation information for genes of interest (72h)
gene_annotations_72h <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype", "go_id", "name_1006", "definition_1006"),
  filters = "hgnc_symbol",
  values = similar_genes_72h, #sostituisci "gene_column" con il nome della colonna dei geni
  mart = ensembl
)

# Filter annotations for specific pathways
immune_annotations_72h <- gene_annotations_72h[grepl
                  ("immune", gene_annotations_72h$name_1006, ignore.case = TRUE), 
                 c("hgnc_symbol", "gene_biotype", "name_1006", "definition_1006")]
immune_annotations_72h <- unique(immune_annotations_72h)


lipid_annotations_72h <- gene_annotations_72h[grepl
                  ("lipid", gene_annotations_72h$name_1006, ignore.case = TRUE), 
                 c("hgnc_symbol", "gene_biotype", "name_1006", "definition_1006")]
lipid_annotations_72h <- unique(lipid_annotations_72h)


mitochondria_annotations_72h <- gene_annotations_72h[grepl 
            ("mitochondria", gene_annotations_72h$name_1006, ignore.case = TRUE),
           c("hgnc_symbol", "gene_biotype", "name_1006", "definition_1006")]
mitochondria_annotations_72h <- unique(mitochondria_annotations_72h)

# Combine all annotations into a single data set
combined_annotations_72h <- bind_rows(
  immune_annotations_72h,
  lipid_annotations_72h,
  mitochondria_annotations_72h
) %>%
  distinct()  # Remove duplicate rows

# Save to Excel
write.xlsx(combined_annotations_72h, "ensembl_72h.xlsx", rowNames = FALSE)


####Decrease####
decrease_data <- increase_decrease[increase_decrease$Change_24h == "Decrease", ]
#download 
write.xlsx(decrease_data, "decrease_data.xlsx", rowNames=TRUE)
