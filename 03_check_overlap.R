# Load required libraries
library(data.table)
library(plyr)

# Load all datasets into a named list
datasets <- list(
  a = fread("PC_Kings_2016_2017_Psychv1.1_Hg19.bim", header = FALSE),
  b = fread("PI_ESALQ_2020_GSAv2.0_Hg19.bim", header = FALSE),
  c = fread("PSC_Kings_2019_GSAv3.0_Hg19.bim", header = FALSE),
  d = fread("PSC_Kings_2016_2019_Psychv1.1_Hg19.bim", header = FALSE),
  e = fread("TP_USP_2020_GSAv1.0_Hg19.bim", header = FALSE),
  f = fread("SC_CHOP_2017_OmniExpressv1.0_Hg19.bim", header = FALSE),
  g = fread("PS_CHOP_2017_OmniExpressv1.1_Hg19.bim", header = FALSE)
)

# Generate all pairwise combinations of datasets
pairwise_combinations <- combn(names(datasets), 2, simplify = FALSE)

# Perform inner joins for all pairwise combinations
pairwise_results <- lapply(pairwise_combinations, function(pair) {
  join_all(list(datasets[[pair[1]]], datasets[[pair[2]]]), by = "V2", type = "inner")
})

# Name the results based on the dataset pairs
names(pairwise_results) <- sapply(pairwise_combinations, function(pair) paste(pair, collapse = ""))

# Combine individual datasets and pairwise results into a single list
all_results <- c(datasets, pairwise_results)

# Count the number of rows in each dataset and pairwise result
row_counts <- sapply(all_results, nrow)

# Print the row counts
print(row_counts)