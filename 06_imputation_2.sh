# Define variables path and name
plink_file=/path/to/plink_file
liftover=/path/to/liftover
path_file=/path/to/path_file
output=/path/to/output

# Step 1: Exclude problematic regions (novel CUPs) from the dataset and convert chromosomes to PLINK format.
plink --bfile "$plink_file" --exclude range "$liftover"/ALL_GRCh37_novel_CUPs_plink.bed --output-chr M --make-bed --out "$plink_file"_del1

# Step 2: Recode the dataset to standard PLINK format after excluding problematic regions.
plink --bfile "$plink_file"_del1 --recode --output-chr M --out "$plink_file"

# Step 3: Perform genome build conversion from Hg19 to Hg38 using the liftOver tool.
# The liftOverPlink.py script converts the MAP file to the new genome build using the specified chain file.
python3 "$liftover"/liftOverPlink.py --bin "$liftover"/liftOver --map "$plink_file".map --out "$path_file"/lifted --chain "$liftover"/hg19ToHg38.over.chain.gz

# Step 4: Remove bad lifts (positions that could not be mapped correctly) using the rmBadLifts.py script.
# This script generates a new MAP file with only successfully lifted positions and logs the bad lifts.
python3 "$liftover"/rmBadLifts.py --map "$path_file"/lifted.map --out "$path_file"/good_lifted.map --log "$path_file"/bad_lifted.dat

# Step 5: Extract the list of SNPs to exclude based on bad lifts.
# The bad_lifted.dat file contains SNPs that failed the liftOver process.
cut -f 2 "$path_file"/bad_lifted.dat > "$path_file"/to_exclude.dat

# Step 6: Add unlifted SNPs from the BED file to the exclusion list.
# The unlifted SNPs are extracted from the BED file and appended to the exclusion list.
cut -f 4 "$path_file"/lifted.bed.unlifted | sed "/^#/d" >> "$path_file"/to_exclude.dat 

# Step 7: Exclude the bad SNPs from the original dataset and recode it to PLINK format.
plink --file "$plink_file" --recode --out "$path_file"/lifted --exclude "$path_file"/to_exclude.dat 

# Step 8: Create a new dataset using the good lifted MAP file.
# This ensures the final dataset is aligned with the new genome build (Hg38).
plink --ped "$path_file"/lifted.ped --map "$path_file"/good_lifted.map --recode --out "$output"

# Step 9: Convert the final dataset to binary PLINK format (BED/BIM/FAM) for efficient storage and analysis.
plink --ped "$output".ped --map "$output".map --make-bed --out "$output"