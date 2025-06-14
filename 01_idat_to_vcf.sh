# Add the user's bin directory to the PATH environment variable for executing custom scripts or binaries.
export PATH="$HOME/bin:$PATH"

# Set the BCFTOOLS_PLUGINS environment variable to point to the user's bin directory, where bcftools plugins are located.
export BCFTOOLS_PLUGINS="$HOME/bin"

# Define variables path and name
path_to_idat_folder=/path/to/idat_folder
bpm_manifest_file=/path/to/bpm_manifest_file
egt_cluster_file=/path/to/egt_cluster_file
path_to_gtc_folder=/path/to/gtc_folder
csv_manifest_file=/path/to/csv_manifest_file
ref=/path/to/reference_genome
out_prefix=/path/to/output_prefix
sample_sheet_file=/path/to/sample_sheet_file
sex_file=/path/to/sex_file

# Check information for IDAT files (e.g., date, chip version, etc.)
# This command uses the `bcftools +gtc2vcf` plugin to inspect IDAT files.
bcftools +gtc2vcf \
  -i -g "$path_to_idat_folder"

# Step 1: Convert IDAT files to GTC format.
# The `bcftools +idat2gtc` plugin is used for this conversion.
# Required inputs include the BPM manifest file, EGT cluster file, and the folder containing IDAT files.
# The output is written to the specified GTC folder.
bcftools +idat2gtc \
  --bpm "$bpm_manifest_file" \
  --egt "$egt_cluster_file" \
  --idats "$path_to_idat_folder" \
  --output "$path_to_gtc_folder"

# Step 2: Convert GTC files to VCF format.
# The `bcftools +gtc2vcf` plugin is used for this conversion.
# Required inputs include the BPM manifest file, CSV manifest file, EGT cluster file, GTC folder, and a reference genome.
# Additional outputs include a TSV file with extra information and a VCF file.
bcftools +gtc2vcf \
  --bpm "$bpm_manifest_file" \
  --csv "$csv_manifest_file" \
  --egt "$egt_cluster_file" \
  --gtcs "$path_to_gtc_folder" \
  --fasta-ref "$ref" \
  --extra "$out_prefix.tsv" | \

# Step 3: Reheader the VCF file to update sample names using a sample sheet.
  bcftools reheader --samples "$sample_sheet_file" | \

# Step 4: Sort the VCF file to ensure proper ordering of variants.
  bcftools sort -Ou -T ./bcftools. | \

# Step 5: Normalize the VCF file to ensure consistent representation of variants.
# This includes left-aligning indels and splitting multiallelic sites.
  bcftools norm --threads 15 --no-version -Oz -c x -f "$ref" | \

# Step 6: Save the final VCF file and index it for efficient querying.
  tee "$out_prefix.vcf.gz" | \
  bcftools index --threads 15 -ft --output "$out_prefix.vcf.gz.tbi"

# Step 7: Convert the VCF file to PLINK format (bed, bim, fam files).
# PLINK is a widely used tool for genetic data analysis.
# The `--const-fid` option ensures constant family IDs, and `--update-sex` updates sex information from a file.
plink --vcf "$out_prefix.vcf.gz" --const-fid --update-sex "$sex_file" --make-bed --out "$out_prefix"