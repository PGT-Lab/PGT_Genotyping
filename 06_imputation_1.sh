# Define variables path and name
data=/path/to/data
data_freq=/path/to/data_freq
check=/path/to/check
output_path=/path/to/output_path
prep_TOPMed=/path/to/prep_TOPMed
data_path=/path/to/data_path
summary_file=/path/to/summary_file
batch=/path/to/batch

# Step 1: Calculate allele frequencies for the dataset and save the output.
plink --freq --bfile "$data" --out "$data_freq"

# Step 2: Run the HRC-1000G-check-bim.pl script to check and prepare the BIM file for imputation.
# This script checks the dataset against the reference panel and identifies issues like allele mismatches.
perl "$check"/HRC-1000G-check-bim.pl -b "$data".bim -f "$data"_freq.frq -r "$check"/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz -h

# Step 3: Execute the generated PLINK commands to fix the dataset.
bash Run-plink.sh

# Step 4: Create output directories for the processed data and preparation files.
mkdir -p "$output_path" "$prep_TOPMed"

# Step 5: Handle the X chromosome by renaming files and replacing "chr23" with "chrX".
find "$data_path" -type f -name "*chr23*" | while IFS= read -r file; do 
  new_name="${file//chr23/chrX}"  # Replace "chr23" with "chrX" in the filename.
  mv "$file" "$new_name"         # Rename the file.
done

# Step 6: Convert VCF files for chromosomes 1-22 and X to PLINK-compatible format.
for i in {1..22} X; do 
  plink --vcf "$data-updated-chr${i}.vcf" --output-chr chrMT --recode vcf --out "$data-updated-chr${i}" --noweb
done

# Step 7: Check and fix reference mismatches in the VCF files using bcftools +fixref.
# This step ensures the dataset matches the reference genome, flipping alleles if necessary.
for i in {1..22} X; do 
  bcftools +fixref "$data"-updated-chr"$i".vcf -Oz -o "$data"-updated-chr"$i".ref.vcf -- \
    -d -f "$check"/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -m flip 2>&1 | tee "$data"-updated-chr"$i".ref.vcf.log
done 

# Step 8: Summarize unresolved variants for each chromosome and calculate the total unresolved variants.
echo -e "Chromosome\tUnresolved_Variants" > "$data_path"/"$summary_file"
total_unresolved=0
for i in {1..22} X; do
  log_file="$data"-updated-chr"$i".ref.vcf.log
  if [ -f "$log_file" ]; then
    # Extract the count of unresolved variants from the log file.
    unresolved_count=$(grep -oP "^NS\s+unresolved\s+\K\d+" "$log_file")
    unresolved_count=${unresolved_count:-0}  # Default to 0 if no unresolved variants are found.
    total_unresolved=$((total_unresolved + unresolved_count))  # Update the total unresolved count.
    echo -e chr"$i"\t"$unresolved_count" >> "$data_path/"$summary_file"
  else
    echo -e chr"$i\tN/A" >> "$data_path/"$summary_file"  # Log "N/A" if the log file is missing.
  fi
done
echo -e Total\t"$total_unresolved" >> "$data_path"/"$summary_file"  # Add the total unresolved count to the summary.
cat "$data_path"/"$summary_file"  # Display the summary file.

# Step 9: Sort and compress the fixed VCF files for each chromosome.
for i in {1..22} X; do 
  vcf-sort "$data"-updated-chr"$i".ref.vcf | bgzip -c > "$output_path"/"$batch"_chr"$i".vcf.gz
done

# Step 10: Extract files from zipped archives for each chromosome using 7z.
for i in {1..22} X; do 
  7z e "$output_path"/chr_"$i".zip -p'password' -o"$output_path"/unzip
done