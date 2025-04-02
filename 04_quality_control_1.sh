# Extract phenotype information (columns 1, 2, and 6) and save it to a file for updating PLINK files.
awk '{print $1, $2, $6}' "$fenotype_information" > "$feno_uptade_file".txt

# Merge phenotype information with the genotyped PLINK dataset.
plink --bfile "$plink_files_genotyped" --pheno "$feno_uptade_file".txt --pheno-merge --make-bed --out "$plink_files_genotyped"_pheno

# Update individual IDs (columns 1, 2, and 3) in the phenotype file and save it for updating PLINK files.
awk '{print $1, $2, $3, $2}' "$fenotype_information" > "$feno_uptade_file".txt

# Update family and individual IDs in the PLINK dataset.
plink --bfile "$plink_files_genotyped"_pheno --update-ids "$feno_uptade_file".txt --make-bed --out "$plink_files_genotyped"_pheno_FID

# Perform initial quality control (QC) by filtering SNPs with:
# - Minor allele frequency (MAF) < 0.01
# - Genotyping rate < 98% (geno 0.02)
# - Hardy-Weinberg equilibrium (HWE) p-value < 1e-10
plink --bfile "$plink_files_genotyped" --maf 0.01 --geno 0.02 --hwe 1e-10 --make-bed --out "$plink_files_genotyped"_QC

# Filter individuals with missing genotype rates > 2% (mind 0.02).
plink --bfile "$plink_files_genotyped"_QC --mind 0.02 --make-bed --out "$plink_files_genotyped"_QC_mind

# Keep only SNPs with valid A, C, G, T alleles (no ambiguous SNPs).
plink --bfile "$plink_files_genotyped"_QC_mind --snps-only just-acgt --make-bed --out "$plink_files_genotyped"_QC_mind_ACTG

# Calculate identity-by-descent (IBD) statistics for relatedness analysis, excluding pairs with IBD > 0.2.
plink --bfile "$plink_files_genotyped"_QC_mind_ACTG --genome --min 0.2 --out "$plink_files_genotyped"_ibd

# Exclude sex chromosomes (XY) from the dataset for further analysis.
plink --bfile "$plink_files_genotyped"_ibd --not-chr XY --make-bed --out "$plink_files_genotyped"_semXY

# Calculate missing genotype rates for individuals.
plink --bfile "$plink_files_genotyped"_semXY --missing --out indiv_missing

# Calculate heterozygosity rates for individuals.
plink --bfile "$plink_files_genotyped"_semXY --het --out hetero

# Run the R script to identify individuals with:
# - High missing genotype rates (> 3%)
# - Extreme heterozygosity (> 3 standard deviations from the mean)
Rscript 04_quality_control_2.R

# Remove individuals identified as outliers in the R script (fail-imisshet-qc.txt).
plink --bfile "$plink_files_genotyped"_semXY --remove fail-imisshet-qc.txt --make-bed --out "$plink_files_genotyped"_QCed