###############################
### 1. Load Required Libraries
###############################
# (Uncomment installation lines if packages are not installed)
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "org.Hs.eg.db", "DESeq2", "apeglm", "pheatmap", "vsn"))
# install.packages(c("dplyr", "RColorBrewer", "readr", "openxlsx", "factoextra"))
# Install the package if not already installed
install.packages("sf")

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
library(readxl)


setwd("C:\\Users\\hp\\Desktop\\QuPath\\Analisi in R")


#####################################################
###### 6 months #####################################
#####################################################

# Define the file path
df1 <- read.csv("Data_MO35_6m.csv", sep=";", encoding="ISO-8859-1")
df2 <- read.csv("Data_MO34_6m.csv", sep=";", encoding="ISO-8859-1")
df3 <- read.csv("Data_MO33_6m.csv", sep=";", encoding="ISO-8859-1")

head(df1)  # Show the first 6 rows
str(df1)   # Show the structure of the columns
summary(df1)  # Statistics on the columns

df1 <- read.csv("Data_MO35_6m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df1[, 4:ncol(df1)] <- lapply(df1[, 4:ncol(df1)], as.numeric)  # Convert numeric columns
summary(df1)
df1[, 4:ncol(df1)] <- lapply(df1[, 4:ncol(df1)], function(x) as.numeric(gsub(",", ".", x)))
df1 <- na.omit(df1)

df2 <- read.csv("Data_MO34_6m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df2[, 4:ncol(df2)] <- lapply(df2[, 4:ncol(df2)], as.numeric)  # Convert numeric columns
summary(df2)
df2[, 4:ncol(df2)] <- lapply(df2[, 4:ncol(df2)], function(x) as.numeric(gsub(",", ".", x)))
df2 <- na.omit(df2)

df3 <- read.csv("Data_MO33_6m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df3[, 4:ncol(df3)] <- lapply(df3[, 4:ncol(df3)], as.numeric)  # Convert numeric columns
summary(df3)
df3[, 4:ncol(df3)] <- lapply(df3[, 4:ncol(df3)], function(x) as.numeric(gsub(",", ".", x)))
df3 <- na.omit(df3)

colnames(df2) <- colnames(df1)
colnames(df3) <- colnames(df1)

#################################
# Combine the datasets
df_tot <- rbind(df1, df2, df3)

# Remove NA values
df_6months<- na.omit(df_tot)

# Save data frame as CSV
write.csv(df_6months, "df_6months.csv")

# Calculate the average Nucleus..Area and Cell..Area for each image and each classification
summary_data_6months <- df_6months %>%
  group_by(Image, Classification) %>%
  summarise(
    Mean_Nucleus_Area = mean(Nucleus..Area, na.rm = TRUE),
    Mean_Cell_Area = mean(Cell..Area, na.rm = TRUE),
    Count = n()
  )

# Calculate the average of the averages and the total count for each classification
final_summary_6months <- summary_data_6months %>%
  group_by(Classification) %>%
  summarise(
    Mean_Nucleus_Area_Total = mean(Mean_Nucleus_Area, na.rm = TRUE),  # Media delle medie di Nucleus..Area
    Mean_Cell_Area_Total = mean(Mean_Cell_Area, na.rm = TRUE),  # Media delle medie di Cell..Area
    Mean_Count_Total = mean(Count, na.rm = TRUE)  # Media di Count
  )

# Function to compute averages for each dataset
calculate_summary <- function(df, MO_name) {
  summary_data <- df %>%
    group_by(Image, Classification) %>%
    summarise(
      Mean_Nucleus_Area = mean(Nucleus..Area, na.rm = TRUE),
      Mean_Cell_Area = mean(Cell..Area, na.rm = TRUE),
      Count = n()
    )
  
  final_summary <- summary_data %>%
    group_by(Classification) %>%
    summarise(
      Mean_Nucleus_Area_Total = mean(Mean_Nucleus_Area, na.rm = TRUE),
      Mean_Cell_Area_Total = mean(Mean_Cell_Area, na.rm = TRUE),
      Mean_Count_Total = mean(Count, na.rm = TRUE)
    ) %>%
    mutate(MO = MO_name)  # Add a column indicating the dataset source
  
  return(final_summary)
}

# Compute summaries for each MO
summary_MO35 <- calculate_summary(df1, "MO35")
summary_MO34 <- calculate_summary(df2, "MO34")
summary_MO33 <- calculate_summary(df3, "MO33")

# Combine all summaries into one final dataset
final_summary_all_MO <- bind_rows(summary_MO35, summary_MO34, summary_MO33)


#####################################################
###### wild type ####################################
#####################################################

# Define the file path
df_wt1<- read.csv("Data_MO25_WT.csv", sep=";", encoding="ISO-8859-1")
df_wt2<- read.csv("Data_MO26_WT.csv", sep=";", encoding="ISO-8859-1")
df_wt3<- read.csv("Data_MO27_WT.csv", sep=";", encoding="ISO-8859-1")
df_wt4<- read.csv("Data_MO28_WT.csv", sep=";", encoding="ISO-8859-1")


df_wt1 <- read.csv("Data_MO25_WT.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df_wt1[, 4:ncol(df_wt1)] <- lapply(df_wt1[, 4:ncol(df_wt1)], as.numeric)  # Convert numeric columns
summary(df_wt1)
df_wt1[, 4:ncol(df_wt1)] <- lapply(df_wt1[, 4:ncol(df_wt1)], function(x) as.numeric(gsub(",", ".", x)))
df_wt1 <- na.omit(df_wt1)

df_wt2 <- read.csv("Data_MO26_WT.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df_wt2[, 4:ncol(df_wt2)] <- lapply(df_wt2[, 4:ncol(df_wt2)], as.numeric)  # Convert numeric columns
summary(df_wt2)
df_wt2[, 4:ncol(df_wt2)] <- lapply(df_wt2[, 4:ncol(df_wt2)], function(x) as.numeric(gsub(",", ".", x)))
df_wt2 <- na.omit(df_wt2)

df_wt3 <- read.csv("Data_MO27_WT.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df_wt3[, 4:ncol(df_wt3)] <- lapply(df_wt3[, 4:ncol(df_wt3)], as.numeric)  # Convert numeric columns
summary(df_wt3)
df_wt3[, 4:ncol(df_wt3)] <- lapply(df_wt3[, 4:ncol(df_wt3)], function(x) as.numeric(gsub(",", ".", x)))
df_wt3 <- na.omit(df_wt3)

df_wt4 <- read.csv("Data_MO28_WT.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df_wt4[, 4:ncol(df_wt4)] <- lapply(df_wt4[, 4:ncol(df_wt4)], as.numeric)  # Convert numeric columns
summary(df_wt4)
df_wt4[, 4:ncol(df_wt4)] <- lapply(df_wt4[, 4:ncol(df_wt4)], function(x) as.numeric(gsub(",", ".", x)))
df_wt4 <- na.omit(df_wt4)

##################################################
# Combine the wild-type data
df_wt <- rbind(df_wt1, df_wt2, df_wt3, df_wt4)

# Make sure the column names match for the wild-type data (same as the other datasets)
colnames(df_wt) <- colnames(df1)

# Save data frame as CSV
write.csv(df_wt, "df_wt.csv")

# Group by 'Image' and 'Classification' to calculate the mean Nucleus Area, Cell::Area, and mean Count per image (for Iba1)
summary_wt <- df_wt %>%
  filter(Classification == "Iba1") %>%  # Ensure we only work with Iba1 cells
  group_by(Image, Classification) %>%   # Keep Classification for final summary
  summarise(
    Mean_Nucleus_Area = mean(Nucleus..Area, na.rm = TRUE),
    Mean_Cell_Area = mean(Cell..Area, na.rm = TRUE),  # Add mean Cell..Area
    Mean_Count = n()  # Count the number of occurrences (cells) for each image
  )

# Now, calculate the average of the means (Mean_Nucleus_Area, Mean_Cell_Area, and Mean_Count) across all images, keeping Classification
final_summary_wt <- summary_wt %>%
  group_by(Classification) %>%
  summarise(
    Mean_Nucleus_Area_Total = mean(Mean_Nucleus_Area, na.rm = TRUE),
    Mean_Cell_Area_Total = mean(Mean_Cell_Area, na.rm = TRUE),  # Calculate mean of Mean_Cell_Area
    Mean_Count_Total = mean(Mean_Count, na.rm = TRUE)
  )

# Add Amylo-Glow row with NA values for Mean_Nucleus_Area, Mean_Cell_Area, and Mean_Count (assuming no values for Amylo-Glow)
amylo_glow_row <- tibble(
  Classification = "Amylo-Glow",
  Mean_Nucleus_Area_Total = NA_real_,
  Mean_Cell_Area_Total = NA_real_,  # NA for Cell::Area
  Mean_Count_Total = 0  # Assuming there are no counts for Amylo-Glow
)

# Bind this row with the final summary
final_summary_wt <- bind_rows(final_summary_wt, amylo_glow_row)

# final_summary_wt now contains both wild-type data and Amylo-Glow, with the average of Cell..Area included

# Compute summaries for each WT dataset
summary_MO25 <- calculate_summary(df_wt1, "MO25")
summary_MO26 <- calculate_summary(df_wt2, "MO26")
summary_MO27 <- calculate_summary(df_wt3, "MO27")
summary_MO28 <- calculate_summary(df_wt4, "MO28")

# Combine all summaries into one final dataset
final_summary_all_WT <- bind_rows(summary_MO25, summary_MO26, summary_MO27, summary_MO28)


#####################################################
###### 4 months #####################################
#####################################################


# Define the file path
df4 <- read.csv("Data_MO20_4m.csv", sep=";", encoding="ISO-8859-1")
df5 <- read.csv("Data_MO31_4m.csv", sep=";", encoding="ISO-8859-1")
df6 <- read.csv("Data_MO32_4m.csv", sep=";", encoding="ISO-8859-1")

df4 <- read.csv("Data_MO20_4m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df4[, 4:ncol(df4)] <- lapply(df4[, 4:ncol(df4)], as.numeric)  # Convert numeric columns
summary(df4)
df4[, 4:ncol(df4)] <- lapply(df4[, 4:ncol(df4)], function(x) as.numeric(gsub(",", ".", x)))
df4 <- na.omit(df4)

df5 <- read.csv("Data_MO31_4m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df5[, 4:ncol(df5)] <- lapply(df5[, 4:ncol(df5)], as.numeric)  # Convert numeric columns
summary(df5)
df5[, 4:ncol(df5)] <- lapply(df5[, 4:ncol(df5)], function(x) as.numeric(gsub(",", ".", x)))
df5 <- na.omit(df5)

df6 <- read.csv("Data_MO32_4m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df6[, 4:ncol(df6)] <- lapply(df6[, 4:ncol(df6)], as.numeric)  # Convert numeric columns
summary(df6)
df6[, 4:ncol(df6)] <- lapply(df6[, 4:ncol(df6)], function(x) as.numeric(gsub(",", ".", x)))
df6 <- na.omit(df6)

#############################################################################
# Check column names for each dataset
colnames(df4) <- colnames(df1)
colnames(df5) <- colnames(df4)

# Combine the 4-month datasets
df_4_months <- rbind(df4, df5, df6)

# Remove NA values
df_4_months <- na.omit(df_4_months)

# Save data frame as CSV
write.csv(df_4_months, "df_4_months.csv")


#Calculate the average Nucleus..Area for each image and each classification
summary_data_4months <- df_4_months %>%
  group_by(Image, Classification) %>%
  summarise(
    Mean_Nucleus_Area = mean(Nucleus..Area, na.rm = TRUE),
    Mean_Cell_Area = mean(Cell..Area, na.rm = TRUE),
    Count = n()
  )

#Calculate the average of the averages and the total count for each classification
final_summary_4months <- summary_data_4months %>%
  group_by(Classification) %>%
  summarise(
    Mean_Nucleus_Area_Total = mean(Mean_Nucleus_Area, na.rm = TRUE),  # Media delle medie di Nucleus..Area
    Mean_Cell_Area_Total = mean(Mean_Cell_Area, na.rm = TRUE),  # Media delle medie di Cell..Area
    Mean_Count_Total = mean(Count, na.rm = TRUE)  # Media di Count
  )

# Compute summaries for each dataset
summary_MO20 <- calculate_summary(df4, "MO20")
summary_MO31 <- calculate_summary(df5, "MO31")
summary_MO32 <- calculate_summary(df6, "MO32")

# Combine all summaries into one final dataset
final_summary_all_4m <- bind_rows(summary_MO20, summary_MO31, summary_MO32)

#####################################################
###### 2 months #####################################
#####################################################


# Define the file path
df7 <- read.csv("Data_MO18_2m.csv", sep=";", encoding="ISO-8859-1")
df8 <- read.csv("Data_MO29_2m.csv", sep=";", encoding="ISO-8859-1")
df9 <- read.csv("Data_MO30_2m.csv", sep=";", encoding="ISO-8859-1")

df7 <- read.csv("Data_MO18_2m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df7[, 4:ncol(df7)] <- lapply(df7[, 4:ncol(df7)], as.numeric)  # Convert numeric columns
summary(df7)
df7[, 4:ncol(df7)] <- lapply(df7[, 4:ncol(df7)], function(x) as.numeric(gsub(",", ".", x)))
df7 <- na.omit(df7)

df8 <- read.csv("Data_MO29_2m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df8[, 4:ncol(df8)] <- lapply(df8[, 4:ncol(df8)], as.numeric)  # Convert numeric columns
summary(df8)
df8[, 4:ncol(df8)] <- lapply(df8[, 4:ncol(df8)], function(x) as.numeric(gsub(",", ".", x)))
df8 <- na.omit(df8)

df9 <- read.csv("Data_MO30_2m.csv", sep=";", encoding="ISO-8859-1", stringsAsFactors=FALSE)
df9[, 4:ncol(df9)] <- lapply(df9[, 4:ncol(df9)], as.numeric)  # Convert numeric columns
summary(df9)
df9[, 4:ncol(df9)] <- lapply(df9[, 4:ncol(df9)], function(x) as.numeric(gsub(",", ".", x)))
df9 <- na.omit(df9)

#############################################################################
# Check column names for each dataset
colnames(df7) <- colnames(df1)
colnames(df8) <- colnames(df7)
colnames(df9) <- colnames(df7)

# Combine the 2-month datasets
df_2_months <- rbind(df7, df8, df9)

# Remove NA values
df_2_months <- na.omit(df_2_months)

# Save data frame as CSV
write.csv(df_2_months, "df_2_months.csv")

# Group by 'Image' and 'Classification' to calculate the mean Nucleus Area, Cell::Area, and mean Count per image (for Iba1)
summary_2months <- df_2_months %>%
  filter(Classification == "Iba1") %>%  # Ensure we only work with Iba1 cells
  group_by(Image, Classification) %>%   # Keep Classification for final summary
  summarise(
    Mean_Nucleus_Area = mean(Nucleus..Area, na.rm = TRUE),
    Mean_Cell_Area = mean(Cell..Area, na.rm = TRUE),  # Add mean Cell..Area
    Mean_Count = n()  # Count the number of occurrences (cells) for each image
  )

# Now, calculate the average of the means (Mean_Nucleus_Area, Mean_Cell_Area, and Mean_Count) across all images, keeping Classification
final_summary_2months <- summary_2months %>%
  group_by(Classification) %>%
  summarise(
    Mean_Nucleus_Area_Total = mean(Mean_Nucleus_Area, na.rm = TRUE),
    Mean_Cell_Area_Total = mean(Mean_Cell_Area, na.rm = TRUE),  # Calculate mean of Mean_Cell_Area
    Mean_Count_Total = mean(Mean_Count, na.rm = TRUE)
  )

# Add Amylo-Glow row with NA values for Mean_Nucleus_Area, Mean_Cell_Area, and Mean_Count (assuming no values for Amylo-Glow)
amylo_glow_row_2months <- tibble(
  Classification = "Amylo-Glow",
  Mean_Nucleus_Area_Total = NA_real_,
  Mean_Cell_Area_Total = NA_real_,  # NA for Cell..Area
  Mean_Count_Total = 0  # Assuming there are no counts for Amylo-Glow
)

# Bind this row with the final summary
final_summary_2months <- bind_rows(final_summary_2months, amylo_glow_row_2months)

# final_summary_2months now contains both wild-type data and Amylo-Glow, with the average of Cell..Area included


# Compute summaries for each dataset
summary_MO18 <- calculate_summary(df7, "MO18")
summary_MO29 <- calculate_summary(df8, "MO29")
summary_MO30 <- calculate_summary(df9, "MO30")

# Combine all summaries into one final dataset
final_summary_all_2m <- bind_rows(summary_MO18, summary_MO29, summary_MO30)

#############################################################################
########### ggplot ##########################################################
#############################################################################

# Rename 'Count' to 'Mean_Count' in 4M and 6M datasets
colnames(summary_data_4months)[colnames(summary_data_4months) == "Count"] <- "Mean_Count"
colnames(summary_data_6months)[colnames(summary_data_6months) == "Count"] <- "Mean_Count"

# Add the 'group' column to each dataset
summary_wt$group <- "WT"
summary_2months$group <- "2M"
summary_data_4months$group <- "4M"
summary_data_6months$group <- "6M"

# Combine the datasets
combined_data <- rbind(summary_wt, summary_2months, summary_data_4months, summary_data_6months)

# Save data frame as CSV
write.csv(combined_data, "Combined_data.csv")

# View the combined dataframe
head(combined_data)

# Filter the data to include only rows where Classification == "Iba1"
iba1_data <- combined_data[combined_data$Classification == "Iba1", ]

# Convert 'group' to a factor with specific order
iba1_data$group <- factor(iba1_data$group, levels = c("WT", "2M", "4M", "6M"))

# Create the plot
ggplot(iba1_data, aes(x = group, y = Mean_Nucleus_Area)) +
  # Bar plot (with means for each group)
  stat_summary(fun = "mean", geom = "bar", fill = "skyblue", color = "skyblue", width = 0.6, position = "dodge") +
  # Scatter plot (individual data points)
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  # Add lines (whiskers) from the min to the max for each group
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  # Customize the plot
  labs(x = "Group", y = "Mean Nucleus Area", title = "Mean Nucleus Area by Group for Iba1") +
  theme_minimal() +
  theme(legend.position = "none")
############################################################################################

# Create the plot
ggplot(iba1_data, aes(x = group, y = Mean_Cell_Area)) +
  # Bar plot (with means for each group)
  stat_summary(fun = "mean", geom = "bar", fill = "lightgreen", color = "lightgreen", width = 0.6, position = "dodge") +
  # Scatter plot (individual data points)
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  # Add lines (whiskers) from the min to the max for each group
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  # Customize the plot
  labs(x = "Group", y = "Mean Microglia size", title = "Mean Microglia size") +
  theme_minimal() +
  theme(legend.position = "none")
############################################################################################
# Create the plot
ggplot(iba1_data, aes(x = group, y = Mean_Count)) +
  # Bar plot (with means for each group)
  stat_summary(fun = "mean", geom = "bar", fill = "violet", color = "violet", width = 0.6, position = "dodge") +
  # Scatter plot (individual data points)
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  # Add lines (whiskers) from the min to the max for each group
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  # Customize the plot
  labs(x = "Group", y = "Mean Count", title = "Microglia Count") +
  theme_minimal() +
  theme(legend.position = "none")

############################################################################################
# Filter the data to include only rows where Classification == "Amylo-Glow"
amylo_glow_data <- combined_data[combined_data$Classification == "Amylo-Glow", ]

# Add empty rows for missing groups to ensure they appear on the x-axis
missing_data <- data.frame(
  group = c("WT", "2M"),
  Mean_Nucleus_Area = NA,  # Missing values for Mean_Nucleus_Area
  Classification = "Amylo-Glow"
)

# Combine the missing data with the original data
combined_amylo_glow_data <- rbind(amylo_glow_data, missing_data)

# Convert 'group' to a factor with specific order
combined_amylo_glow_data$group <- factor(combined_amylo_glow_data$group, levels = c("WT", "2M", "4M", "6M"))

# Create the plot
ggplot(combined_amylo_glow_data, aes(x = group, y = Mean_Cell_Area)) +
  # Bar plot (with means for each group)
  stat_summary(fun = "mean", geom = "bar", fill = "lightcoral", color = "lightcoral", width = 0.6, position = "dodge") +
  # Scatter plot (individual data points) with color mapped to group
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  # Add lines (whiskers) from the min to the max for each group
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  # Customize the plot
  labs(x = "Group", y = "Mean Plaque Area", title = "Plaque Size") +
  theme_minimal() +
  theme(legend.position = "none")

###########################################################################################

# Add missing groups for Mean_Count
missing_data_count <- data.frame(
  group = c("WT", "2M"),
  Mean_Count = NA,  # Missing values for Mean_Count
  Classification = "Amylo-Glow"
)

# Combine with the original Amylo-Glow data
combined_amylo_glow_count <- rbind(amylo_glow_data, missing_data_count)

# Convert 'group' to a factor to enforce order
combined_amylo_glow_count$group <- factor(combined_amylo_glow_count$group, levels = c("WT", "2M", "4M", "6M"))

# Create the plot for Mean_Count
ggplot(combined_amylo_glow_count, aes(x = group, y = Mean_Count)) +
  # Bar plot (means for each group)
  stat_summary(fun = "mean", geom = "bar", fill = "lightcoral", color = "lightcoral", width = 0.6, position = "dodge") +
  # Scatter plot (individual data points)
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  # Error bars (whiskers from min to max)
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  # Customize the plot
  labs(x = "Group", y = "Mean Plaque Count", title = "Plaque Count") +
  theme_minimal() +
  theme(legend.position = "none")

# Load necessary packages
library(ggplot2)
library(patchwork)

# Store each plot in a variable

# 1. Mean Nucleus Area for Iba1
plot_iba1_nucleus <- ggplot(iba1_data, aes(x = group, y = Mean_Nucleus_Area)) +
  stat_summary(fun = "mean", geom = "bar", fill = "skyblue", color = "skyblue", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Nucleus Area", title = "Mean Nucleus Area") +
  theme_minimal() + theme(legend.position = "none")

# 2. Mean Cell Area for Iba1
plot_iba1_cell <- ggplot(iba1_data, aes(x = group, y = Mean_Cell_Area)) +
  stat_summary(fun = "mean", geom = "bar", fill = "lightgreen", color = "lightgreen", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Cell Area", title = "Mean Cell Area") +
  theme_minimal() + theme(legend.position = "none")

# 3. Mean Count for Iba1
plot_iba1_count <- ggplot(iba1_data, aes(x = group, y = Mean_Count)) +
  stat_summary(fun = "mean", geom = "bar", fill = "violet", color = "violet", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Count", title = "Microglia Count") +
  theme_minimal() + theme(legend.position = "none")

# 4. Mean Cell Area for Amylo-Glow
plot_amylo_glow_cell <- ggplot(combined_amylo_glow_data, aes(x = group, y = Mean_Cell_Area)) +
  stat_summary(fun = "mean", geom = "bar", fill = "lightcoral", color = "lightcoral", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Plaque Area", title = "Amylo-Glow: Plaque Size") +
  theme_minimal() + theme(legend.position = "none")

# 5. Mean Count for Amylo-Glow
plot_amylo_glow_count <- ggplot(combined_amylo_glow_count, aes(x = group, y = Mean_Count)) +
  stat_summary(fun = "mean", geom = "bar", fill = "orange", color = "orange", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Plaque Count", title = "Amylo-Glow: Plaque Count") +
  theme_minimal() + theme(legend.position = "none")

# Arrange the plots in a 2x3 layout
final_plot <- (plot_iba1_nucleus | plot_iba1_cell | plot_iba1_count) /
  (plot_amylo_glow_cell | plot_amylo_glow_count)

# Save the final combined plot as an image
ggsave("all_plots_combined.png", final_plot, width = 14, height = 10, dpi = 300)

###########################################################################
########## BAR PLOT WITH MOUSE ###################################
###########################################################################
library(ggrepel)

# Add the 'group' column to each dataset
final_summary_all_WT$group <- "WT"
final_summary_all_2m$group <- "2M"
final_summary_all_4m$group <- "4M"
final_summary_all_MO$group <- "6M"

# Combine the datasets
combined_data2 <- rbind(final_summary_all_WT, final_summary_all_2m, final_summary_all_4m, final_summary_all_MO)

# View the combined dataframe
head(combined_data2)

# Ensure 'group' is a factor
combined_data2$group <- factor(combined_data2$group, levels = c("WT", "2M", "4M", "6M"))


# Filter the data to include only rows where Classification == "Iba1"
iba1_data2 <- combined_data2[combined_data2$Classification == "Iba1", ]


# Create the plot for Mean_Nucleus_Area_Total - iba1
ggplot(iba1_data2, aes(x = group, y = Mean_Nucleus_Area_Total)) +
  # Bar plot with means for each group
  stat_summary(fun = "mean", geom = "bar", fill = "skyblue", color = "skyblue", width = 0.6) +
  # Scatter plot (individual data points)
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("WT" = "red", "2M" = "blue", "4M" = "darkgreen", "6M" = "black")) +
  # Add the dots with text for MO values
  geom_text(aes(label = MO), size = 3, vjust = -0.5, color = "black") +
  # Error bars
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Nucleus Area", title = "Microglia") +
  theme_minimal() +
  theme(legend.position = "none")

# Create the plot for Mean_Cell_Area_Total - iba1
ggplot(iba1_data2, aes(x = group, y = Mean_Cell_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "lightgreen", color = "lightgreen", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("WT" = "red", "2M" = "blue", "4M" = "darkgreen", "6M" = "black")) +
  geom_text(aes(label = MO), size = 3, vjust = -0.5, color = "black") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Cell Area", title = "Microglia") +
  theme_minimal() +
  theme(legend.position = "none")

# Create the plot for Mean_Count_Total
ggplot(iba1_data2, aes(x = group, y = Mean_Count_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "violet", color = "violet", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("WT" = "red", "2M" = "blue", "4M" = "darkgreen", "6M" = "black")) +
  geom_text(aes(label = MO), size = 3, vjust = -0.5, color = "black") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Count", title = "Microglia Count by Group") +
  theme_minimal() +
  theme(legend.position = "none")

#####################################################################################

# Filter the data to include only rows where Classification == "Amylo-Glow"
amylo_glow_data2 <- combined_data2[combined_data2$Classification == "Amylo-Glow", ]

# Create missing data for WT and 2M groups with the same column structure as amylo_glow_data2
missing_data_amylo_glow <- data.frame(
  Classification = rep("Amylo-Glow", 2),  # Set Classification as "Amylo-Glow"
  Mean_Nucleus_Area_Total = rep(NA, 2),   # Set NA for Mean_Nucleus_Area_Total
  Mean_Cell_Area_Total = rep(NA, 2),      # Set NA for Mean_Cell_Area_Total
  Mean_Count_Total = rep(NA, 2),          # Set NA for Mean_Count_Total
  MO = rep(NA, 2),                        # Set NA for MO
  group = factor(c("WT", "2M"), levels = c("WT", "2M", "4M", "6M"))  # Specify the groups (WT, 2M)
)

# Ensure the missing data has the same column names and structure as amylo_glow_data2
missing_data_amylo_glow <- missing_data_amylo_glow[, names(amylo_glow_data2)]

# Combine the original data with the missing data
combined_amylo_glow_data2 <- rbind(amylo_glow_data2, missing_data_amylo_glow)

# Now let's check if the combined data frame looks correct
str(combined_amylo_glow_data2)

# Create the plot for Mean_Cell_Area_Total - Amylo-Glow
ggplot(combined_amylo_glow_data2, aes(x = group, y = Mean_Cell_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "lightgreen", color = "lightgreen", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("WT" = "black", "2M" = "darkgreen", "4M" = "blue", "6M" = "red")) +
  geom_text(aes(label = MO), size = 3, vjust = -0.5, color = "black") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Plaque Size", title = "Plaque Size (Including WT and 2M)") +
  theme_minimal() +
  theme(legend.position = "none")

# Create the plot for Mean_Count_Total - Amylo-Glow
ggplot(combined_amylo_glow_data2, aes(x = group, y = Mean_Count_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "violet", color = "violet", width = 0.6) +
  geom_jitter(aes(color = group), size = 2, width = 0.2, alpha = 0.6) +
  scale_color_manual(values = c("WT" = "black", "2M" = "darkgreen", "4M" = "blue", "6M" = "red")) +
  geom_text(aes(label = MO), size = 3, vjust = -0.5, color = "black") +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Count", title = "Plaque Count (Including WT and 2M)") +
  theme_minimal() +
  theme(legend.position = "none")

########################################################################
# save plots 
# Iba1 plots
plot_iba1_nucleus_area <- ggplot(iba1_data2, aes(x = group, y = Mean_Nucleus_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "skyblue", color = "skyblue", width = 0.6) +
  geom_jitter(aes(color = group), size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +
  scale_color_manual(values = c("WT" = "red", "2M" = "blue", "4M" = "darkgreen", "6M" = "black")) +
  ggrepel::geom_text_repel(aes(label = MO), size = 3, color = "black", max.overlaps = 10) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Nucleus Area", title = "Microglia Nucleus Area") +
  theme_minimal() + theme(legend.position = "none")

plot_iba1_cell_area <- ggplot(iba1_data2, aes(x = group, y = Mean_Cell_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "lightgreen", color = "lightgreen", width = 0.6) +
  geom_jitter(aes(color = group), size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +
  scale_color_manual(values = c("WT" = "red", "2M" = "blue", "4M" = "darkgreen", "6M" = "black")) +
  ggrepel::geom_text_repel(aes(label = MO), size = 3, color = "black", max.overlaps = 10) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Cell Area", title = "Microglia Cell Area") +
  theme_minimal() + theme(legend.position = "none")

plot_iba1_count <- ggplot(iba1_data2, aes(x = group, y = Mean_Count_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "violet", color = "violet", width = 0.6) +
  geom_jitter(aes(color = group), size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +
  scale_color_manual(values = c("WT" = "red", "2M" = "blue", "4M" = "darkgreen", "6M" = "black")) +
  ggrepel::geom_text_repel(aes(label = MO), size = 3, color = "black", max.overlaps = 10) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Count", title = "Microglia Count by Group") +
  theme_minimal() + theme(legend.position = "none")

# Amylo-Glow plots
plot_amylo_glow_cell_area <- ggplot(combined_amylo_glow_data2, aes(x = group, y = Mean_Cell_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "lightgreen", color = "lightgreen", width = 0.6) +
  geom_jitter(aes(color = group), size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +
  scale_color_manual(values = c("WT" = "black", "2M" = "darkgreen", "4M" = "blue", "6M" = "red")) +
  ggrepel::geom_text_repel(aes(label = MO), size = 3, color = "black", max.overlaps = 10) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Plaque Size", title = "Plaque Size (Including WT and 2M)") +
  theme_minimal() + theme(legend.position = "none")

plot_amylo_glow_count <- ggplot(combined_amylo_glow_data2, aes(x = group, y = Mean_Count_Total)) +
  stat_summary(fun = "mean", geom = "bar", fill = "violet", color = "violet", width = 0.6) +
  geom_jitter(aes(color = group), size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +
  scale_color_manual(values = c("WT" = "black", "2M" = "darkgreen", "4M" = "blue", "6M" = "red")) +
  ggrepel::geom_text_repel(aes(label = MO), size = 3, color = "black", max.overlaps = 10) + 
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  labs(x = "Group", y = "Mean Count", title = "Plaque Count (Including WT and 2M)") +
  theme_minimal() + theme(legend.position = "none")

# Load the gridExtra package
library(gridExtra)

# Combine the plots in a 3x2 grid (adjust as necessary)
combined_plots <- grid.arrange(
  plot_iba1_nucleus_area, plot_iba1_cell_area,
  plot_iba1_count, plot_amylo_glow_cell_area,
  plot_amylo_glow_count, ncol = 2
)

# Save the combined plots to a PNG file
ggsave("combined_bar_plots.png", combined_plots, width = 12, height = 10, dpi = 600)



#####################################################################################
################# ANOVA #############################################################
#####################################################################################

library(ggsignif)

# ANOVA for Iba1 data (Microglia)
anova_iba1_nucleus <- aov(Mean_Nucleus_Area_Total ~ group, data = iba1_data2)
anova_iba1_cell_area <- aov(Mean_Cell_Area_Total ~ group, data = iba1_data2)
anova_iba1_count <- aov(Mean_Count_Total ~ group, data = iba1_data2)

# ANOVA for Amylo-Glow data (Amyloid Plaques)
anova_amylo_glow_cell_area <- aov(Mean_Cell_Area_Total ~ group, data = combined_amylo_glow_data2)
anova_amylo_glow_count <- aov(Mean_Count_Total ~ group, data = combined_amylo_glow_data2)

# Tukey HSD for pairwise comparisons (if ANOVA is significant)
tukey_iba1_nucleus <- TukeyHSD(anova_iba1_nucleus)
tukey_iba1_cell_area <- TukeyHSD(anova_iba1_cell_area)
tukey_iba1_count <- TukeyHSD(anova_iba1_count)

tukey_amylo_glow_cell_area <- TukeyHSD(anova_amylo_glow_cell_area)
tukey_amylo_glow_count <- TukeyHSD(anova_amylo_glow_count)

# For the Iba1 Nucleus Area plot:
plot_iba1_nucleus <- ggplot(iba1_data2, aes(x = group, y = Mean_Nucleus_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", aes(fill = group, color = group), width = 0.6) +  # Bars and borders colored by group
  geom_jitter(color = "black", size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +  # Dots all black
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, color = "black") +  # Error bars
  geom_signif(comparisons = list(c("2M", "4M"), c("2M", "6M"), c("4M", "6M")), map_signif_level = TRUE) +  # Significance annotations
  labs(x = "Group", y = "Mean Nucleus Area", title = "Microglia Nucleus Area") +
  scale_fill_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Custom fill colors for bars
  scale_color_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Matching borders with fill
  theme_minimal() + theme(legend.position = "none")

# For the Iba1 Cell Area plot:
plot_iba1_cell_area <- ggplot(iba1_data2, aes(x = group, y = Mean_Cell_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", aes(fill = group, color = group), width = 0.6) +  # Bars and borders colored by group
  geom_jitter(color = "black", size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +  # Dots all black
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, color = "black") +  # Error bars
  geom_signif(comparisons = list(c("2M", "4M"), c("2M", "6M"), c("4M", "6M")), map_signif_level = TRUE) +  # Significance annotations
  labs(x = "Group", y = "Mean Cell Area", title = "Microglia Cell Area") +
  scale_fill_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Custom fill colors for bars
  scale_color_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Matching borders with fill
  theme_minimal() + theme(legend.position = "none")

# For the Iba1 Count plot:
plot_iba1_count <- ggplot(iba1_data2, aes(x = group, y = Mean_Count_Total)) +
  stat_summary(fun = "mean", geom = "bar", aes(fill = group, color = group), width = 0.6) +  # Bars and borders colored by group
  geom_jitter(color = "black", size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +  # Dots all black
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, color = "black") +  # Error bars
  geom_signif(comparisons = list(c("2M", "4M"), c("2M", "6M"), c("4M", "6M")), map_signif_level = TRUE) +  # Significance annotations
  labs(x = "Group", y = "Mean Count", title = "Microglia Count by Group") +
  scale_fill_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Custom fill colors for bars
  scale_color_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Matching borders with fill
  theme_minimal() + theme(legend.position = "none")


# Amylo-Glow plots
# Plot without MOxx label, with age color for bars, and black dots, using mean_se for error bars

plot_amylo_glow_cell_area <- ggplot(combined_amylo_glow_data2, aes(x = group, y = Mean_Cell_Area_Total)) +
  stat_summary(fun = "mean", geom = "bar", aes(fill = group, color = group), width = 0.6) +  # Bars and borders colored by group
  geom_jitter(color = "black", size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +  # Dots all black
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, color = "black") +  # Error bars
  geom_signif(comparisons = list(c("2M", "4M"), c("2M", "6M"), c("4M", "6M")), map_signif_level = TRUE) +  # Significance annotations
  labs(x = "Group", y = "Mean Plaque Size", title = "Plaque Size (Including WT and 2M)") +
  scale_fill_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Custom fill colors for bars
  scale_color_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Matching borders with fill
  theme_minimal() + theme(legend.position = "none")

plot_amylo_glow_count <- ggplot(combined_amylo_glow_data2, aes(x = group, y = Mean_Count_Total)) +
  stat_summary(fun = "mean", geom = "bar", aes(fill = group, color = group), width = 0.6) +  # Bars and borders colored by group
  geom_jitter(color = "black", size = 2.5, width = 0.2, alpha = 0.6, stroke = 0.3) +  # Dots all black
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2, color = "black") +  # Error bars
  geom_signif(comparisons = list(c("2M", "4M"), c("2M", "6M"), c("4M", "6M")), map_signif_level = TRUE) +  # Significance annotations
  labs(x = "Group", y = "Mean Count", title = "Plaque Count (Including WT and 2M)") +
  scale_fill_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Custom fill colors for bars
  scale_color_manual(values = c("WT" = "#A7C7E7", "2M" = "pink", "4M" = "lightgreen", "6M" = "#D1B2FF")) +  # Matching borders with fill
  theme_minimal() + theme(legend.position = "none")



# Combine the plots in a 3x2 grid (adjust as necessary)
combined_plots_all <- grid.arrange(
  plot_iba1_nucleus, plot_iba1_cell_area,
  plot_iba1_count, plot_amylo_glow_cell_area,
  plot_amylo_glow_count, ncol = 2
)


# Save the combined plots to a PNG file
ggsave("combined_bar_plots.png", combined_plots_all, width = 12, height = 10, dpi = 600)
