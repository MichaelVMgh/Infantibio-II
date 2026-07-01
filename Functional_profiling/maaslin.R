#BiocManager::install("Maaslin2") #To install Maaslin2
library(Maaslin2)
library(tidyverse)
library(vegan)
library(data.table)

# set Working Directory (WDir)
setwd("C:/path/to/your/data/analyses/")
# Load the Humann3 unstratified abundances data table (gene pathways or gene families)
df_input_data = read.delim("genepathways_unstratified_2batches_babies_without_24mo_transposed.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
# df_input_data_genus = read.table("merged_abundance_table_genus_transposed.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
# df_input_data_babies = df_input_data[df_input_data$TimePoint != "mother",]
# # show the first 5 rows and columns
# # each row is a sample, each column is a species
# df_input_data[1:5, 1:5]
# df_input_data_3mo = read.table("merged_abundance_table_species_transposed_3mo.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
# 
# Load the metadata
df_input_metadata = read.delim("mapping_file_combined_babies_without_24mo.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
# df_input_metadata[1:5, ]
# df_input_metadata_3mo = read.table("mapping_file_3mo.txt", header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
# df_input_metadata_babies = df_input_metadata[df_input_metadata["TimePoint"] != "mom",]



# Subset df_input_data to keep only the rows (samples) that are also in df_input_metadata
df_input_data <- df_input_data[rownames(df_input_data) %in% rownames(df_input_metadata), ]


# Assign factors
df_input_metadata$Feeding <- factor(df_input_metadata$Feeding, levels=c("BF", "formula"))
#df_input_metadata$HistoryAtTimePoint <- factor(df_input_metadata$HistoryAtTimePoint, levels=c("BF", "A", "B", "C", "D", "infant_24mo"))
df_input_metadata$HistoryAtTimePoint <- factor(df_input_metadata$HistoryAtTimePoint, levels=c("A", "B", "C", "D"))
#df_input_metadata$Clade <- factor(df_input_metadata$Clade, levels=c("1", "2"))
#df_input_metadata$HistoryAtTimePoint_combined <- factor(df_input_metadata$HistoryAtTimePoint_combined, levels=c("breastfeeding", "A", "B", "C", "D"))
df_input_metadata$TypeOfBirth <- factor(df_input_metadata$TypeOfBirth, levels=c("vaginal", "c-section"))
#df_input_metadata$AntibioticsMother <- factor(df_input_metadata$AntibioticsMother, levels=c("yes", "no"))
#df_input_metadata$Sex <- factor(df_input_metadata$Sex, levels=c("male", "female"))
df_input_metadata$TimePoint <- factor(df_input_metadata$TimePoint, levels=c("1mo", "3mo", "7mo", "12mo"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex, levels=c("male", "female"))
df_input_metadata$AntibioticsMother <- factor(df_input_metadata$AntibioticsMother, levels=c("yes", "no"))
df_input_metadata$Supplementation <- factor(df_input_metadata$Supplementation, levels=c("Yes", "No"))

max_value <- max(unlist(df_input_data), na.rm = TRUE)
#summary(df_input_data)
output_dir = "C:/path/to/your/data/analyses/"

#Run Maaslin
fit_data_1 = Maaslin2(
  input_data = df_input_data,
  input_metadata = df_input_metadata,
  output = output_dir,
  fixed_effects = c("Feeding","TypeOfBirth","Sex","AntibioticsMother"), #Choose covariates
  reference = c("Feeding,BF","TypeOfBirth,vaginal","Sex,female","AntibioticsMother,no"),
  random_effects = c("ID"),
  min_abundance = 1000,
  normalization = "NONE",
  transform = "LOG",
  analysis_method = "LM",
  max_significance = 0.25,
  correction = "BH")





