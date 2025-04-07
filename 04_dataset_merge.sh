# Define variables path and name
gsa1=/path/to/gsa1
gsa2=/path/to/gsa2
gsa3=/path/to/gsa3
gsa4=/path/to/gsa4
gsa5=/path/to/gsa5
gsa1_new=/path/to/gsa1_new
gsa2_new=/path/to/gsa2_new
gsa3_new=/path/to/gsa3_new
gsa4_new=/path/to/gsa4_new
gsa5_new=/path/to/gsa5_new
gsa_path=/path/to/gsa_path
gsa_merged=/path/to/gsa_merged

# Step 1: Remove duplicate variants from each dataset and generate a list of SNPs
plink --bfile "$gsa1" --list-duplicate-vars --out "$gsa1_new"  # Identify duplicate variants in dataset 1
plink --bfile "$gsa1" --exclude "$gsa1_new".dupvar --make-bed --out "$gsa1_new"_nodup  # Exclude duplicate variants
plink --bfile "$gsa1_new"_nodup --write-snplist --out "$gsa1_new"_nodup  # Write the list of SNPs
sort -o "$gsa1_new"_nodup.snplist "$gsa1_new"_nodup.snplist  # Sort the SNP list

# Repeat the same process for datasets 2 to 5
plink --bfile "$gsa2" --list-duplicate-vars --out "$gsa2_new"
plink --bfile "$gsa2" --exclude "$gsa2_new".dupvar --make-bed --out "$gsa2_new"_nodup
plink --bfile "$gsa2_new"_nodup --write-snplist --out "$gsa2_new"_nodup
sort -o "$gsa2_new"_nodup.snplist "$gsa2_new"_nodup.snplist

plink --bfile "$gsa3" --list-duplicate-vars --out "$gsa3_new"
plink --bfile "$gsa3" --exclude "$gsa3_new".dupvar --make-bed --out "$gsa3_new"_nodup
plink --bfile "$gsa3_new"_nodup --write-snplist --out "$gsa3_new"_nodup
sort -o "$gsa3_new"_nodup.snplist "$gsa3_new"_nodup.snplist

plink --bfile "$gsa4" --list-duplicate-vars --out "$gsa4_new"
plink --bfile "$gsa4" --exclude "$gsa4_new".dupvar --make-bed --out "$gsa4_new"_nodup
plink --bfile "$gsa4_new"_nodup --write-snplist --out "$gsa4_new"_nodup
sort -o "$gsa4_new"_nodup.snplist "$gsa4_new"_nodup.snplist

plink --bfile "$gsa5" --list-duplicate-vars --out "$gsa5_new"
plink --bfile "$gsa5" --exclude "$gsa5_new".dupvar --make-bed --out "$gsa5_new"_nodup
plink --bfile "$gsa5_new"_nodup --write-snplist --out "$gsa5_new"_nodup
sort -o "$gsa5_new"_nodup.snplist "$gsa5_new"_nodup.snplist

# Step 2: Identify common SNPs across all datasets
comm -12 "$gsa1_new"_nodup.snplist "$gsa2_new"_nodup.snplist | \
comm -12 - "$gsa3_new"_nodup.snplist | \
comm -12 - "$gsa4_new"_nodup.snplist | \
comm -12 - "$gsa5_new"_nodup.snplist > "$gsa_merged".snplist  # Find common SNPs across all datasets

# Step 3: Filter each dataset to include only the common SNPs
plink --bfile "$gsa1_new"_nodup --extract "$gsa_merged".snplist --make-bed --out "$gsa1_new"_nodup_common
plink --bfile "$gsa2_new"_nodup --extract "$gsa_merged".snplist --make-bed --out "$gsa2_new"_nodup_common
plink --bfile "$gsa3_new"_nodup --extract "$gsa_merged".snplist --make-bed --out "$gsa3_new"_nodup_common
plink --bfile "$gsa4_new"_nodup --extract "$gsa_merged".snplist --make-bed --out "$gsa4_new"_nodup_common
plink --bfile "$gsa5_new"_nodup --extract "$gsa_merged".snplist --make-bed --out "$gsa5_new"_nodup_common

# Step 4: Update FID (Family IDs) for datasets 2 to 5
awk '{print $1, $2, "1", $2}' "$gsa2_new"_nodup_common.fam > "$gsa2_new"_updatefids.txt  # Update FID for dataset 2
plink --bfile "$gsa2_new"_nodup_common --update-ids "$gsa2_new"_updatefids.txt --make-bed --out "$gsa2_new"_nodup_common_fid1

awk '{print $1, $2, "2", $2}' "$gsa3_new"_nodup_common.fam > "$gsa3_new"_updatefids.txt  # Update FID for dataset 3
plink --bfile "$gsa3_new"_nodup_common --update-ids "$gsa3_new"_updatefids.txt --make-bed --out "$gsa3_new"_nodup_common_fid2

awk '{print $1, $2, "3", $2}' "$gsa4_new"_nodup_common.fam > "$gsa4_new"_updatefids.txt  # Update FID for dataset 4
plink --bfile "$gsa4_new"_nodup_common --update-ids "$gsa4_new"_updatefids.txt --make-bed --out "$gsa4_new"_nodup_common_fid3

awk '{print $1, $2, "4", $2}' "$gsa5_new"_nodup_common.fam > "$gsa5_new"_updatefids.txt  # Update FID for dataset 5
plink --bfile "$gsa5_new"_nodup_common --update-ids "$gsa5_new"_updatefids.txt --make-bed --out "$gsa5_new"_nodup_common_fid4

# Step 5: Merge all datasets sequentially
plink --bfile "$gsa1_new"_nodup_common --bmerge "$gsa2_new"_nodup_common_fid1 --make-bed --out "$gsa_path"/temp1  # Merge dataset 1 and 2
plink --bfile "$gsa_path"/temp1 --bmerge "$gsa3_new"_nodup_common_fid2 --make-bed --out "$gsa_path"/temp2  # Merge with dataset 3
plink --bfile "$gsa_path"/temp2 --bmerge "$gsa4_new"_nodup_common_fid3 --make-bed --out "$gsa_path"/temp3  # Merge with dataset 4
plink --bfile "$gsa_path"/temp3 --bmerge "$gsa5_new"_nodup_common_fid4 --make-bed --out "$gsa_merged"  # Merge with dataset 5