# Load required data
imiss <- read.table('indiv_missing.imiss', header = TRUE)
het <- read.table('hetero.het', header = TRUE)

# Calculate the proportion of heterozygosity
het$P_HET <- (het$N.NM. - het$O.HOM.) / het$N.NM.

# Define thresholds for heterozygosity based on 3 standard deviations
upper_3sd <- mean(het$P_HET) + 3 * sd(het$P_HET)
lower_3sd <- mean(het$P_HET) - 3 * sd(het$P_HET)

# Generate a plot of missing genotypes vs heterozygosity
pdf('imiss-vs-het.pdf')
plot(
  log10(imiss$F_MISS), het$P_HET,
  xlab = 'log10(Proportion of missing genotypes)',
  ylab = 'Proportion Heterozygous',
  xlim = c(-4, 0), ylim = c(0, 0.5)
)
axis(side = 1, labels = FALSE)
mtext(c(-4, -3, -2, -1, 0), side = 1, at = c(-4, -3, -2, -1, 0), line = 1)
abline(h = upper_3sd, col = 'red', lty = 2)  # Upper threshold
abline(h = lower_3sd, col = 'red', lty = 2)  # Lower threshold
abline(v = log10(0.03), col = 'red', lty = 2)  # Missingness threshold
dev.off()

# Identify individuals to remove based on thresholds
imiss_rem <- subset(imiss, imiss$F_MISS > 0.03)[, 1:2]  # Individuals with high missingness
het_rem <- subset(het, het$P_HET > upper_3sd | het$P_HET < lower_3sd)[, 1:2]  # Outliers in heterozygosity

# Combine individuals to remove into a single list
indiv_rem <- rbind(imiss_rem, het_rem)

# Save the list of individuals to remove
write.table(
  indiv_rem, 'fail-imisshet-qc.txt',
  col.names = FALSE, row.names = FALSE, quote = FALSE, sep = '\t'
)