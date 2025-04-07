# Define variables path and name
plink_files_genotyped=/path/to/plink_files_genotyped

# Perform initial quality control (QC) by filtering SNPs with:
# - Minor allele frequency (MAF) < 0.01
# - Genotyping rate < 98% (geno 0.02)
# - Hardy-Weinberg equilibrium (HWE) p-value < 1e-10
plink --bfile "$plink_files_genotyped" --maf 0.01 --geno 0.02 --hwe 1e-10 --make-bed --out "$plink_files_genotyped"_QC

# Filter individuals with missing genotype rates > 2% (mind 0.02).
plink --bfile "$plink_files_genotyped"_QC --mind 0.02 --make-bed --out "$plink_files_genotyped"_QC_mind

# Keep only SNPs with valid A, C, G, T alleles.
plink --bfile "$plink_files_genotyped"_QC_mind --snps-only just-acgt --make-bed --out "$plink_files_genotyped"_QC_mind_ACTG

# Calculate identity-by-descent (IBD) statistics for relatedness analysis, excluding pairs with IBD > 0.2.
plink --bfile "$plink_files_genotyped"_QC_mind_ACTG --genome --min 0.2 --out "$plink_files_genotyped"_ibd

# Exclude sex chromosomes (XY) from the dataset for further analysis.
plink --bfile "$plink_files_genotyped"_QC_mind_ACTG --not-chr XY --make-bed --out "$plink_files_genotyped"_semXY

# Calculate missing genotype rates for individuals.
plink --bfile "$plink_files_genotyped"_semXY --missing --out indiv_missing

# Calculate heterozygosity rates for individuals.
plink --bfile "$plink_files_genotyped"_semXY --het --out hetero

# Run the R script to identify individuals with:
# - High missing genotype rates (> 3%)
# - Extreme heterozygosity (> 3 standard deviations from the mean)
Rscript 04_quality_control_2.R

# Remove individuals identified as outliers in the R script (fail-imisshet-qc.txt).
plink --bfile "$plink_files_genotyped"_ibd --remove fail-imisshet-qc.txt --make-bed --out "$plink_files_genotyped"_semhet

# Perform a sex check to identify discrepancies between reported and genetic sex.
plink --bfile "$plink_files_genotyped"_semhet --sex-check --out check_XY

# Extract individuals with sex discrepancies (marked as "PROBLEM") from the sex check output.
grep PROBLEM check_XY.sexcheck | awk '{print $1, $2}' > check_sex_fail.txt

# Remove individuals with sex discrepancies from the dataset and create the final quality-controlled dataset.
plink --bfile "$plink_files_genotyped"_semhet --remove check_sex_fail.txt --make-bed --out "$plink_files_genotyped"_QCed