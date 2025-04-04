# Define the path to the directory containing the unzipped imputed VCF files
output_path_GSA=/path/to/unziped/imputed/file/vcf/

# Find all .dose.vcf.gz files in the specified directory, sort them numerically, 
# and filter variants with R2 > 0.8 using bcftools. The filtered files are saved with a new name.
find "${output_path_GSA}" -name "*.dose.vcf.gz" | sort -V  | parallel -j 5 'bcftools view -i "R2>.8" -Oz -o {.}.filtered_r2_08.vcf.gz {}'

# Find all filtered VCF files, sort them numerically, and concatenate them into a single VCF file.
find "${output_path_GSA}" -name "*.filtered_r2_08.vcf.gz" | sort -V | xargs bcftools concat -Oz -o "${output_path_GSA}/imputed_GSA_r2_08.vcf.gz"

# Index the concatenated VCF file to enable fast querying
bcftools index -t "${output_path_GSA}/imputed_GSA_r2_08.vcf.gz"

# Convert the concatenated VCF file into PLINK binary format (.bed, .bim, .fam) for downstream analysis
plink --vcf ${output_path_GSA}/imputed_GSA_r2_08.vcf.gz --make-bed --out ${output_path_GSA}/imputed_GSA_r2_08