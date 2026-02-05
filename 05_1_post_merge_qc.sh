#!/bin/bash
# Pre-Imputation QC Pipeline (PLINK-based)
# Usage: ./preimputation_qc.sh <INPUT_PREFIX> <OUTPUT_PREFIX>

# ------------------------------------------------------------------------------
# Initialize
# ------------------------------------------------------------------------------
if [ "$#" -ne 2 ]; then
  echo "Error: Usage: $0 <INPUT_PREFIX> <OUTPUT_PREFIX>"
  exit 1
fi

INPUT=$1
OUTPUT=$2
mkdir -p logs

# ------------------------------------------------------------------------------
# Step 1: Initial IBD Check (Relatedness)
# ------------------------------------------------------------------------------
echo "STEP 1: Calculating initial IBD/PI_HAT..."
plink --bfile $INPUT --genome --min 0.2 --out logs/${INPUT}_initial_ibd 2> logs/step1_ibd.log
sort -rk 10 logs/${INPUT}_initial_ibd.genome > logs/${INPUT}_initial_ibd_sorted.genome

# ------------------------------------------------------------------------------
# Step 2: SNP QC Filters
# ------------------------------------------------------------------------------
echo "STEP 2: Filtering indels and missing data..."
plink --bfile $INPUT \
  --snps-only just-acgt \
  --make-bed \
  --out ${INPUT}_snps_only 2> logs/step2_snps.log

plink --bfile ${INPUT}_snps_only \
  --geno 0.05 \
  --make-bed \
  --out ${INPUT}_geno_filtered 2> logs/step2_geno.log

plink --bfile ${INPUT}_geno_filtered \
  --mind 0.02 \
  --make-bed \
  --out ${INPUT}_mind_filtered 2> logs/step2_mind.log

# ------------------------------------------------------------------------------
# Step 3: Heterozygosity/Missingness Outliers (Calls R Script)
# ------------------------------------------------------------------------------
echo "STEP 3: Detecting heterozygosity/missingness outliers..."
plink --bfile ${INPUT}_mind_filtered \
  --missing \
  --out logs/${INPUT}_indiv_missing 2> logs/step3_missing.log

plink --bfile ${INPUT}_mind_filtered \
  --het \
  --out logs/${INPUT}_het 2> logs/step3_het.log

# Call R script for analysis
Rscript scripts/05_2_heterozygosity_outliers.R ${INPUT} 2> logs/step3_r_het.log

# Remove outliers if file exists
if [ -f "logs/${INPUT}_fail_het_imiss.txt" ]; then
  plink --bfile ${INPUT}_mind_filtered \
    --remove logs/${INPUT}_fail_het_imiss.txt \
    --make-bed \
    --out ${INPUT}_het_filtered 2> logs/step3_remove_het.log
else
  echo "No heterozygosity outliers detected. Proceeding..."
  cp ${INPUT}_mind_filtered.{bed,bim,fam} ${INPUT}_het_filtered.{bed,bim,fam}
fi

# ------------------------------------------------------------------------------
# Step 4: Sex Discrepancy Check (Calls R Script)
# ------------------------------------------------------------------------------
echo "STEP 4: Checking sex discrepancies..."
plink --bfile ${INPUT}_het_filtered \ 
  --split-x b37 \ #ATENTION! IF YOUR DATASET FOR SOME REASON IS ALREADY SPLIT, YOU SHOULD SKIP THIS STEP AND CHANGE THE NAME OF THE FILE IN THE NEXT STEP.
  --make-bed \
  --out ${INPUT}_splitx 2> logs/step4_splitx.log

plink --bfile ${INPUT}_splitx \ #ATENTION! IF YOUR DATASET FOR SOME REASON WAS ALREADY SPLIT, YOU SHOULD CHANGE THE INPUT NAME HERE TO ${INPUT}_het_filtered
  --check-sex \
  --out logs/${INPUT}_sexcheck 2> logs/step4_sexcheck.log

# Call R script for analysis
Rscript scripts/05_3_sex_check.R ${INPUT} 2> logs/step4_r_sexcheck.log

# Remove discrepancies if file exists
if [ -f "logs/${INPUT}_sexcheck_problem_samples.txt" ]; then
  plink --bfile ${INPUT}_het_filtered \
    --remove logs/${INPUT}_sexcheck_problem_samples.txt \
    --make-bed \
    --out ${INPUT}_sex_filtered 2> logs/step4_remove_sex.log
else
  echo "No sex discrepancies detected. Proceeding..."
  cp ${INPUT}_het_filtered.{bed,bim,fam} ${INPUT}_sex_filtered.{bed,bim,fam}
fi

# ------------------------------------------------------------------------------
# Step 5: Final QC Filters (MAF, HWE)
# ------------------------------------------------------------------------------
echo "STEP 5: Applying MAF and HWE filters..."
plink --bfile ${INPUT}_sex_filtered \
  --geno 0.02 \
  --make-bed \
  --out ${INPUT}_geno_final 2> logs/step5_geno.log

plink --bfile ${INPUT}_geno_final \
  --maf 0.01 \
  --make-bed \
  --out ${INPUT}_maf_filtered 2> logs/step5_maf.log

plink --bfile ${INPUT}_maf_filtered \
  --hwe 1e-10 \
  --make-bed \
  --out ${OUTPUT}_pre_imputation 2> logs/step5_hwe.log

# ------------------------------------------------------------------------------
# Step 6: Final IBD Check
# ------------------------------------------------------------------------------
echo "STEP 6: Final IBD/PI_HAT calculation..."
plink --bfile ${OUTPUT}_pre_imputation \
  --genome \
  --min 0.2 \
  --out logs/${OUTPUT}_final_ibd 2> logs/step6_ibd.log

sort -rk 10 logs/${OUTPUT}_final_ibd.genome > logs/${OUTPUT}_final_ibd_sorted.genome

# ------------------------------------------------------------------------------
# Completion
# ------------------------------------------------------------------------------
echo "QC pipeline complete! Final output: ${OUTPUT}_pre_imputation.{bed,bim,fam}"
