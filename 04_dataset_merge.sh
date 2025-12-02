# Define variables path and name
dataset_path=/mnt/rosalind_1/Lagc_scz/WCPG/Costa_Rica/3_Imputation
dataset1_path=$dataset_path
dataset2_path=$dataset_path
dataset1=$dataset1_path/merged_ld_cr_QCed
dataset2=$dataset2_path/merged_fams_QCed
dataset1_new=merged_ld_controls_cr_new
dataset2_new=merged_fams_new
dataset_merged=merged_ld_cr_fams

# Step 1: Remove duplicate variants from each dataset and generate a list of SNPs
plink2 --bfile "$dataset1" --set-all-var-ids @:# --make-bed --out "${dataset1_path}/${dataset1_new}"
plink --bfile "${dataset1_path}/${dataset1_new}" --list-duplicate-vars --out "${dataset1_path}/${dataset1_new}"  # Identify duplicate variants in dataset 1
plink --bfile  "${dataset1_path}/${dataset1_new}" --exclude "${dataset1_path}/${dataset1_new}.dupvar" --make-bed --out "${dataset1_path}/${dataset1_new}_nodup"  # Exclude duplicate variants
plink --bfile "${dataset1_path}/${dataset1_new}_nodup" --write-snplist --out "${dataset1_path}/${dataset1_new}_nodup"  # Write the list of SNPs
sort -o "${dataset1_path}/${dataset1_new}_nodup.snplist" "${dataset1_path}/${dataset1_new}_nodup.snplist"  # Sort the SNP list

#plink --bfile "$dataset2" --set-missing-var-ids @:#:$a:$b --make-bed --out "${dataset2_path}/${dataset2_new}_varid"
plink --bfile "$dataset2" --list-duplicate-vars --out  "${dataset2_path}/${dataset2_new}"  # Identify duplicate variants in dataset 1
plink --bfile "${dataset2_path}/${dataset2_new}" --exclude  "${dataset2_path}/${dataset2_new}.dupvar" --make-bed --out "${dataset2_path}/${dataset2_new}_nodup"  # Exclude duplicate variants
plink --bfile "${dataset2_path}/${dataset2_new}_nodup" --write-snplist --out "${dataset2_path}/${dataset2_new}_nodup"  # Write the list of SNPs
sort -o "${dataset2_path}/${dataset2_new}_nodup.snplist" "${dataset2_path}/${dataset2_new}_nodup.snplist"  # Sort the SNP list

# Step 2: Identify common SNPs across all datasets
comm -12 "${dataset1_path}/${dataset1_new}_nodup.snplist" "${dataset2_path}/${dataset2_new}_nodup.snplist" > "${dataset_path}/${dataset_merged}.snplist"  # Find common SNPs across all datasets

# Step 3: Filter each dataset to include only the common SNPs
plink --bfile  "${dataset1_path}/${dataset1_new}_nodup" --extract "${dataset_path}/${dataset_merged}.snplist" --make-bed --out "${dataset1_path}/${dataset1_new}_nodup_common"
plink --bfile  "${dataset2_path}/${dataset2_new}_nodup" --extract "${dataset_path}/${dataset_merged}.snplist" --make-bed --out "${dataset2_path}/${dataset2_new}_nodup_common"

# Step 4: Update FID (Family IDs) to make sure it can handle duplicated IIDs
awk '{print $1, $2, "1", $2}'  "${dataset2_path}/${dataset2_new}_nodup_common.fam" > "${dataset2_path}/${dataset2_new}_updatefids.txt"  # Update FID for dataset 2
plink --bfile  "${dataset2_path}/${dataset2_new}_nodup_common" --update-ids "${dataset2_path}/${dataset2_new}_updatefids.txt" --make-bed --out "${dataset2_path}/${dataset2_new}_nodup_common_fid1"
# Step 5: Merge all datasets sequentially
plink --bfile "${dataset1_path}/${dataset1_new}_nodup_common" --bmerge "${dataset2_path}/${dataset2_new}_nodup_common_fid1" --make-bed --out "$dataset_path"/"$dataset_merged"  # Merge dataset 1 and 2
plink --bfile "${dataset1_path}/${dataset1_new}_nodup_common" --exclude "${dataset_path}/${dataset_merged}-merge.missnp" --make-bed --out "${dataset1_path}/${dataset1_new}_nodup_common_biallelic"
plink --bfile "${dataset2_path}/${dataset2_new}_nodup_common_fid1" --exclude "${dataset_path}/${dataset_merged}-merge.missnp" --make-bed --out  "${dataset2_path}/${dataset2_new}_nodup_common_fid1_biallelic"
plink --bfile "${dataset1_path}/${dataset1_new}_nodup_common_biallelic" --bmerge "${dataset2_path}/${dataset2_new}_nodup_common_fid1_biallelic" --make-bed --out  "${dataset_path}/${dataset_merged}"

