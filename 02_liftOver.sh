# Step 1: Convert PLINK binary files to text format (PED/MAP) for compatibility with liftOver tools.
plink --bfile "$plink_file" --recode --out "$prefix_hg19"

# Step 2: Use the liftOverPlink.py script to lift over the genetic data from Hg19 to Hg38.
# The script requires the liftOver binary, the input MAP file, and the chain file for the conversion.
python3 liftOverPlink.py --bin ./liftOver --map "$prefix_hg19.map" --out lifted --chain hg19ToHg38.over.chain.gz

# Step 3: Remove bad lifts (positions that could not be mapped correctly) using rmBadLifts.py.
# This script generates a new MAP file with only good lifts and logs the bad lifts.
python3 rmBadLifts.py --map lifted.map --out good_lifted.map --log bad_lifted.dat

# Step 4: Extract the list of SNPs to exclude based on bad lifts.
# The bad_lifted.dat file contains SNPs that failed the liftOver process.
cut -f 2 bad_lifted.dat > to_exclude.dat

# Step 5: Add unlifted SNPs from the BED file to the exclusion list.
# The unlifted SNPs are extracted from the BED file and appended to the exclusion list.
cut -f 4 lifted.bed.unlifted | sed "/^#/d" >> to_exclude.dat 

# Step 6: Exclude the bad SNPs from the original PED/MAP files.
# This creates a new PED/MAP file with only the successfully lifted SNPs.
plink --file "$prefix_hg19" --recode --out lifted --exclude to_exclude.dat 

# Step 7: Create a new PED/MAP file using the good lifted MAP file.
# This ensures the final dataset is aligned with the new genome build (Hg38).
plink --ped lifted.ped --map good_lifted.map --recode --out "$output"

# Step 8: Convert the final PED/MAP files back to binary format (BED/BIM/FAM) for efficient storage and analysis.
plink --ped "$output.ped" --map "$output.map" --make-bed --out "$output"