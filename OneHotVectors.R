library(readxl)
library(dplyr)
library(tidyr)

# input file paths
Ben_data_path <- "Cytogenetics_Baseline.tsv"
Skerget_data_path <- "Skerget_etal_Data.xlsx"

# read in data
df_Bens <- read.delim(Ben_data_path, stringsAsFactors = FALSE)
df_Skerget <- read_excel(Skerget_data_path, sheet = "1A_Patient_features")

# get column names that are in both dataframes
common_cols <- intersect(colnames(df_Bens), colnames(df_Skerget))
# drop Hyperdiploid_Call from common columns
common_cols <- common_cols[common_cols != "Hyperdiploid_Call"]

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
write.table(non_matching_samples, file = "non_matching_samples.tsv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t") # nolint