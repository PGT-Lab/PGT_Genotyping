output_path=/path/to/zipped/file/

######format later
# Step 10: Extract files from zipped archives for each chromosome using 7z.
for i in {1..22} X; do 
  7z e "$output_path"/chr_"$i".zip -p'password' -o"$output_path"/unzip
done

# Define the path to the directory containing the unzipped imputed VCF files

# Find all .dose.vcf.gz files in the specified directory, sort them numerically, 
# and filter variants with R2 > 0.8 using bcftools. The filtered files are saved with a new name.
find "$output_path"/unzip -name "*.dose.vcf.gz" | sort -V  | parallel -j 5 'bcftools view -i "R2>.8" -Oz -o {.}.filtered_r2_08.vcf.gz {}'

# Find all filtered VCF files, sort them numerically, and concatenate them into a single VCF file.
find "$output_path"/unzip -name "*.filtered_r2_08.vcf.gz" | sort -V | xargs bcftools concat -Oz -o ""$output_path"/unzip/imputed_r2_08.vcf.gz"

# Index the concatenated VCF file to enable fast querying
bcftools index -t ""$output_path"/unzip/imputed_r2_08.vcf.gz"

# Convert the concatenated VCF file into PLINK binary format (.bed, .bim, .fam) for downstream analysis
plink --vcf "$output_path"/unzip//imputed_r2_08.vcf.gz --make-bed --out "$output_path"/unzip//imputed_r2_08