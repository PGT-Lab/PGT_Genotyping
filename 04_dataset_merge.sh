# Define variables path and name
dataset1=/path/to/dataset1
dataset2=/path/to/dataset2
dataset1_new=/path/to/dataset1_new
dataset2_new=/path/to/dataset2_new
dataset_path=/path/to/dataset_path
dataset_merged=/path/to/dataset_merged

# Step 1: Remove duplicate variants from each dataset and generate a list of SNPs
plink --bfile "$dataset1" --list-duplicate-vars --out "$dataset1_new"  # Identify duplicate variants in dataset 1
plink --bfile "$dataset1" --exclude "$dataset1_new".dupvar --make-bed --out "$dataset1_new"_nodup  # Exclude duplicate variants
plink --bfile "$dataset1_new"_nodup --write-snplist --out "$dataset1_new"_nodup  # Write the list of SNPs
sort -o "$dataset1_new"_nodup.snplist "$dataset1_new"_nodup.snplist  # Sort the SNP list

plink --bfile "$dataset2" --list-duplicate-vars --out "$dataset2_new"  # Identify duplicate variants in dataset 1
plink --bfile "$dataset2" --exclude "$dataset2_new".dupvar --make-bed --out "$dataset2_new"_nodup  # Exclude duplicate variants
plink --bfile "$dataset2_new"_nodup --write-snplist --out "$dataset2_new"_nodup  # Write the list of SNPs
sort -o "$dataset2_new"_nodup.snplist "$dataset2_new"_nodup.snplist  # Sort the SNP list

# Step 2: Identify common SNPs across all datasets
comm -12 "$dataset1_new"_nodup.snplist "$dataset2_new"_nodup.snplist > "$dataset_merged".snplist  # Find common SNPs across all datasets

# Step 3: Filter each dataset to include only the common SNPs
plink --bfile "$dataset1_new"_nodup --extract "$dataset_merged".snplist --make-bed --out "$dataset1_new"_nodup_common
plink --bfile "$dataset2_new"_nodup --extract "$dataset_merged".snplist --make-bed --out "$dataset2_new"_nodup_common

# Step 4: Update FID (Family IDs) to make sure it can handle duplicated IIDs
awk '{print $1, $2, "1", $2}' "$dataset2_new"_nodup_common.fam > "$dataset2_new"_updatefids.txt  # Update FID for dataset 2
plink --bfile "$dataset2_new"_nodup_common --update-ids "$dataset2_new"_updatefids.txt --make-bed --out "$dataset2_new"_nodup_common_fid1

# Step 5: Merge all datasets sequentially
plink --bfile "$dataset1_new"_nodup_common --bmerge "$dataset2_new"_nodup_common_fid1 --make-bed --out "$dataset_path"/"$dataset_merged"  # Merge dataset 1 and 2