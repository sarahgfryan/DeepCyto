library(readxl)
library(dplyr)
library(tidyr)

# input file paths - need to be in working directory
Ben_data_path <- "Cytogenetics_Baseline.csv"
Skerget_data_path <- "Skerget_etal_Data.xlsx"

# read in data
df_Bens <- read.delim(Ben_data_path, stringsAsFactors = FALSE)
df_Skerget <- read_excel(Skerget_data_path, sheet = "1A_Patient_features")

# get column names that are in both dataframes
common_cols <- intersect(colnames(df_Bens), colnames(df_Skerget))
# OPTIONAL: drop Hyperdiploid_Call from common columns (comment out if you want to keep it)
common_cols <- common_cols[common_cols != "Hyperdiploid_Call"]
# OPTIONAL: drop copy number variants from common columns, they have the prefix "Cp_" (comment out if you want to keep them)
#common_cols <- common_cols[!grepl("^Cp_", common_cols)]

# filter both dataframes to only include common columns
df_Bens <- df_Bens %>% select(all_of(common_cols))
df_Skerget <- df_Skerget %>% select(all_of(common_cols))

# save filtered dataframes to new files
write.table(df_Bens, file = "Filtered_Cytogenetics_Baseline.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
write.table(df_Skerget, file = "Filtered_Skerget_etal_Data.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# find samples where cytogenetic call is the same in both dataframes
matching_samples <- df_Bens %>%
  inner_join(df_Skerget, by = common_cols) 
# save matching samples to a text file
write.table(matching_samples, file = "matching_samples.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

#find sample IDs where cytogenetic call is different in both dataframes
non_matching_sampleIDs <- df_Bens %>%
  anti_join(df_Skerget, by = common_cols) %>% 
  distinct(Patient_ID)

#join df_Bens and df_Skerget to get the full information for non-matching samples
non_matching_samples <- non_matching_sampleIDs %>%
  inner_join(df_Bens, by = "Patient_ID") %>%
  inner_join(df_Skerget, by = "Patient_ID", suffix = c("_Bens", "_Skerget"))

# save non-matching samples to a text file
write.table(non_matching_samples, file = "non_matching_samples.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# start with the concordant samples.. count how many cytogenetic abnormalities
colnames(matching_samples) <- gsub("_Call","",colnames(matching_samples))
colnames(matching_samples) <- gsub("Cp_","",colnames(matching_samples))
colnames(matching_samples) <- gsub("_Tx","",colnames(matching_samples))
num_samples <- colSums(matching_samples[,-1]) # exclude Patient ID column from counting

# now do the same for the Ben data, need to ignore NAs
colnames(df_Bens) <- gsub("_Call","",colnames(df_Bens))
colnames(df_Bens) <- gsub("Cp_","",colnames(df_Bens))
colnames(df_Bens) <- gsub("_Tx","",colnames(df_Bens))
num_samples_Bens <- colSums(df_Bens[,-1], na.rm = TRUE) # exclude Patient ID column from counting
print(num_samples_Bens)

# checking validity of WHSC1 data from Bens data
df_WHSC1 <- read.delim("WHSC1.csv", stringsAsFactors = FALSE)
common_cols <- intersect(colnames(df_Bens), colnames(df_WHSC1))
non_matching <- df_Bens %>%
  anti_join(df_WHSC1, by = common_cols) %>% 
  distinct(Patient_ID)
print(non_matching)