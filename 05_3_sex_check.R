#!/usr/bin/env Rscript
# Checks sex discrepancies
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

data <- read.table(paste0("logs/", prefix, "_sexcheck.sexcheck"), header = TRUE)

# Plot
pdf(paste0("logs/", prefix, "_sexcheck.pdf"))
dotchart(data$F, labels = data$IID, xlab = "X-chromosome heterozygosity (F)")
dev.off()

# Identify problematic samples (PROBLEM status and SNPSEX != 0)
problem_samples <- subset(data, STATUS == "PROBLEM" & SNPSEX != 0)[, c(1, 2)]

if (nrow(problem_samples) > 0) {
  write.table(problem_samples, paste0("logs/", prefix, "_sexcheck_problem_samples.txt"), 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("Sex discrepancies detected. Saved to logs/", prefix, "_sexcheck_problem_samples.txt\n", sep = "")
} else {
  cat("No sex discrepancies found.\n")
}
