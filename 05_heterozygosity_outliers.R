#!/usr/bin/env Rscript
# Detects heterozygosity/missingness outliers
args <- commandArgs(trailingOnly = TRUE)
prefix <- args[1]

# Load data
imiss <- read.table(paste0("logs/", prefix, "_indiv_missing.imiss"), header = TRUE)
het <- read.table(paste0("logs/", prefix, "_het.het"), header = TRUE)
het$P_HET <- (het$N.NM. - het$O.HOM.) / het$N.NM.

# Define outliers (3SD)
upper <- mean(het$P_HET) + 3 * sd(het$P_HET)
lower <- mean(het$P_HET) - 3 * sd(het$P_HET)

# Plot
pdf(paste0("logs/", prefix, "_imiss_vs_het.pdf"))
plot(log10(imiss$F_MISS), het$P_HET, 
     xlab = "log10(Missing rate)", ylab = "Heterozygosity rate",
     xlim = c(-4, 0), ylim = c(0, 0.5))
abline(h = c(lower, upper), col = "red", lty = 2)
abline(v = log10(0.03), col = "red", lty = 2)
dev.off()

# Save samples to exclude
fail_het <- subset(het, P_HET > upper | P_HET < lower)[, c(1, 2)]
fail_imiss <- subset(imiss, F_MISS > 0.03)[, c(1, 2)]
fail_qc <- unique(rbind(fail_het, fail_imiss))

if (nrow(fail_qc) > 0) {
  write.table(fail_qc, paste0("logs/", prefix, "_fail_het_imiss.txt"), 
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  cat("Outliers detected. Saved to logs/", prefix, "_fail_het_imiss.txt\n", sep = "")
} else {
  cat("No heterozygosity/missingness outliers found.\n")
}