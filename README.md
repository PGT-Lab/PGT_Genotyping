# PGT quality control pipeline

### Processing of the genotping array datasets at Laboratory of Integrative Neurosciences (LINC) at UNIFESP.

Here we provide the reasoning of parameters and steps taken in the pipeline. The pipeline was designed to be run in a `Linux` environment, and the scripts are written in `bash` and `R`. The pipeline is divided into several steps, each of which is described below. For more detailed information, please refer to the individual scripts.

**Three types of chips were used in the pipeline:**

* Infinium PsychArray v1.1
* Infinium Global Screening Array v1.0/v2.0/v3.0
* HumanOmniExpress-12 v1.1/v1.0

**Major tools in this pipeline:**

* [PLINK v1.9](https://www.cog-genomics.org/plink/)
* [bcftools v1.21](https://github.com/samtools/bcftools)
* [samtools v1.21](https://github.com/samtools/samtools)

For further enquiries, please refer to the authors listed below:

| Authors | Github | Contact |
| ------- | ------ | ------ |
| Lucas Toshio Ito | https://github.com/lcstoshio | lucas.toshio@unifesp.br |
| Adrielle Martins de Oliveira| https://github.com/MartinsAdrielle | adrielle.martins@unifesp.br |
| Gustavo Satoru Kajitani | https://github.com/gkajitani | g.kajitani@unifesp.br |

---

# 1. IDAT to PLINK

![](https://github.com/ccmaues/pgt_images_github/blob/main/01_idat_to_vcf.png?raw=true)

This script convertes IDAT files to GTC format, then to VCF format, and finally to PLINK format for downstream analysis. It uses `bcftools` plugins for file conversions, reheaders and normalizes the VCF file, indexes it for efficient querying, and converts it to PLINK format.

*Requirements for this step:*

- [gtc2vcf](https://github.com/freeseek/gtc2vcf) - Follow installation instructions from their github
- [Chip manifest and cluster files](https://support.illumina.com/array/downloads.html)
- IDAT or GTC files (Illumina)
- Sample sheet and sex file (optional)

*Main commands:*

### **1. IDAT to GTC conversion**

```{bash}
bcftools +idat2gtc --bpm "$bpm_manifest_file" --egt "$egt_cluster_file" --idats "$path_to_idat_folder" --output "$path_to_gtc_folder
```

### **2. GTC to VCF conversion**

```{bash}
bcftools +gtc2vcf \
  --bpm "$bpm_manifest_file" \
  --csv "$csv_manifest_file" \
  --egt "$egt_cluster_file" \
  --gtcs "$path_to_gtc_folder" \
  --fasta-ref "$ref" \
  --extra "$out_prefix.tsv" | \
  bcftools reheader --samples "$sample_sheet_file" | \
  bcftools sort -Ou -T ./bcftools. | \
  bcftools norm --threads 15 --no-version -Oz -c x -f "$ref" | \
  tee "$out_prefix.vcf.gz" | \
  bcftools index --threads 15 -ft --output "$out_prefix.vcf.gz.tbi"
```

### **3. VCF to PLINK v1.9 files (bed bim fam)**

```{bash}
plink --vcf "$out_prefix" --const-fid --update-sex "$sex_file" --make-bed --out "$out_prefix"
```

---

# 2. Check overlap

We estimated the SNV overlap between different chips and cohorts providing a summary of overlaps across those datasets.

![SNVs per chip version](https://github.com/ccmaues/pgt_images_github/blob/main/chip_variants.png?raw=true)

The table below shows the number of overlapping SNVs between different datasets. The diagonal cells represent the dataset itself, while the other cells show the number of overlapping SNVs between the corresponding datasets.

| Dataset                  | PC_Kings_2016_2017 | PI_ESALQ_2020 | PSC_Kings_2019 | PSC_Kings_2016_2019 | TP_USP_2020 | SC_CHOP_2017 | PS_CHOP_2017 |
|--------------------------|--------------------|---------------|---------------|---------------------|-------------|-------------|-------------|
| **PC_Kings_2016_2017**   | -                  | 109542        | 108461        | 587111              | 107784      | 258365      | 256052      |
| **PI_ESALQ_2020**        | 109542             | -             | 635892        | 109542              | 613312      | 139140      | 138358      |
| **PSC_Kings_2019**       | 108461             | 635892        | -             | 108461              | 589402      | 137739      | 136991      |
| **PSC_Kings_2016_2019**  | 587111             | 109542        | 108461        | -                   | 107784      | 258365      | 256052      |
| **TP_USP_2020**          | 107784             | 613312        | 589402        | 107784              | -           | 133277      | 132530      |
| **SC_CHOP_2017**         | 258365             | 139140        | 137739        | 258365              | 133277      | -           | 719023      |
| **PS_CHOP_2017**         | 256052             | 138358        | 136991        | 256052              | 132530      | 719023      | -           |

*Requirements to run:*

```{r}
install.package("data.table")
install.package("plyr")
```

---

# 3. Quality control 1

![](https://raw.githubusercontent.com/ccmaues/pgt_images_github/refs/heads/main/05_QC.png)

Quality control (QC) is critical in genetic studies to ensure data reliability by removing poor-quality samples and variants. This QC step was applied to avoid **false associations, population stratification bias, or spurious results**, compromising study validity. The parameters values were the same as in RICOPILI pipeline, with the addition of relatedness filter.

<!---
add refs here for QC
again... No prunning??
--->

*Main Commands:*

### **1. SNV filter:**

```{bash}
plink --bfile "$genotyped_data" --maf 0.01 --geno 0.02 --hwe 1e-10  --make-bed --out "$genotyped_data"_QC
```

### **2. Sample filter:**

```{bash}
plink --bfile "$genotyped_data"_QC --mind 0.02 --make-bed --out "$genotyped_data"_QC_mind
```

### **3. Keep only ATGC SNVs:**

```{bash}
plink --bfile "$genotyped_data"_QC_mind --snps-only just-acgt --make-bed --out "$genotyped_data"_QC_mind_ACTG
```

### **4. Idendity by Descent**

```{bash}
plink --bfile "$genotyped_data"_QC_mind_ACTG --genome --min 0.2 --out ibd_calculation
```

### **5. Heterozigosity and missingness filter:**
<!---
Wheres the pruning?
--->

```{bash}
plink --bfile "$genotyped_data"_QC_mind_ACTG --missing --out indiv_missing
plink --bfile "$genotyped_data"_QC_mind_ACTG --het --out hetero
```

```{R}
imiss_rem <- subset(imiss, imiss$F_MISS > 0.03)[, 1:2]
het_rem <- subset(het, het$P_HET > upper_3sd | het$P_HET < lower_3sd)[, 1:2] 
```

```{bash}
plink --bfile "$genotyped_data"_QC_mind_ACTG --remove fail-imisshet-qc.txt --make-bed --out "$genotyped_data"_semhet
```

### **7. sex check:**

```{bash}
plink --bfile "$genotyped_data"_semhet --sex-check --out check_XY
grep PROBLEM check_XY.sexcheck | awk '{print $1, $2}' > check_sex_fail.txt
plink --bfile "$genotyped_data"_semhet --remove check_sex_fail.txt --make-bed --out "$genotyped_data"_QCed
```

<!---
put info later
--->
*SNV and sample removal per step:*
| Dataset | MAF | HWE | Geno | Mind | ATGC SNVs | IBD | Heterozygosity | Missingness | Sex-check |
| -------- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| **PC_Kings_2016_2017** | 257,032 | 7,612 | 3,609 | 0 | 0 | 0 | 0 | 0 | 0 |
| **PI_ESALQ_2020** | 145,917 | 15,852 | 71,940 | 6 | 0 | 0 | 0 | 0 | 0 |
| **PSC_Kings_2019** | 15,293 | 6,003 | 298,548 | 6 | 2,441| 0 | 0 | 0 | 0 |
| **PSC_Kings_2016_2019** | 85,573 | 54,891 | 97,262 | 6 | 1,815 | 0 | 0 | 0 | 0 |
| **TP_USP_2020** | 126,363 | 3 | 35,434 | 0 | 95 | 0 | 0 | 0 | 0 |
| **SC_CHOP_2017** | 21,622 | 16,066 | 9,932 | 2 | 0 | 0 | 0 | 0 | 0 |
| **PS_CHOP_2017** | 21,096 | 756 | 15,491 | 0 | 0 | 0 | 0 | 0 | 0 |

---

# 4. Dataset merge

The merge was performed using only the intersection of SNVs between pairs of datasets. This step can be done with two or more datasets, as long as they have sufficient overlap and are based on the same reference genome."

*Main Command:*

### **1. SNV list:**

```{bash}
plink --bfile "$dataset1" --list-duplicate-vars --out "$dataset1_new"
plink --bfile "$dataset1" --exclude "$dataset1_new".dupvar --make-bed --out "$dataset1_new"_nodup
plink --bfile "$dataset1_new"_nodup --write-snplist --out "$dataset1_new"_nodup
sort -o "$dataset1_new"_nodup.snplist "$dataset1_new"_nodup.snplist

plink --bfile "$dataset2" --list-duplicate-vars --out "$dataset2_new"  # Identify duplicate variants in dataset 1
plink --bfile "$dataset2" --exclude "$dataset2_new".dupvar --make-bed --out "$dataset2_new"_nodup  # Exclude duplicate variants
plink --bfile "$dataset2_new"_nodup --write-snplist --out "$dataset2_new"_nodup  # Write the list of SNPs
sort -o "$dataset2_new"_nodup.snplist "$dataset2_new"_nodup.snplist  # Sort the SNP list
```

### **2. Common SNVs between datasets:**

```{bash}
comm -12 "$dataset1_new"_nodup.snplist "$dataset2_new"_nodup.snplist > "$psych_merged".snplist
plink --bfile "$dataset1_new"_nodup --extract "$psych_merged".snplist --make-bed --out "$dataset1_new"_nodup_common
plink --bfile "$dataset2_new"_nodup --extract "$psych_merged".snplist --make-bed --out "$dataset2_new"_nodup_common
```

### **3. FID update:**

```{bash}
awk '{print $1, $2, "1", $2}' "$dataset2_new"_nodup_common.fam > "$dataset2_new"_updatefids.txt
plink --bfile "$dataset2_new"_nodup_common --update-ids "$dataset2_new"_updatefids.txt --make-bed --out "$dataset2_new"_nodup_common_fid1
```

```{bash}
plink --bfile "$dataset1_new"_nodup_common --bmerge "$dataset2_new"_nodup_common_fid1 --make-bed --out "$dataset_path"/"$dataset_merged"
```

---

# 5. Quality Control 2

After merging multiple genotype datasets (e.g., from different batches or studies), perform QC to ensure consistency.

### Pre-Imputation QC Pipeline

A PLINK-based pipeline for quality control (QC) before genotype imputation, adapted from the Ricopili framework.  
**Key steps**: IBD/PI_HAT checks, SNP/indel filtering, missingness, heterozygosity, sex discrepancy, and HWE/MAF filters.

### Requirements
- PLINK 1.9+
- R (for plots)
- Unix shell (`sort`, `bgzip`, etc.)

### Usage
Run scripts preImputation_QC.sh, heterozygosity_outliers.R and sex_check.R.

### **1. Initial IBD/PI_HAT check (relatedness):**
```bash
plink --bfile $INPUT --genome --min 0.2 --out ${INPUT}_initial_ibd
sort -rk 10 ${INPUT}_initial_ibd.genome > ${INPUT}_initial_ibd_sorted.genome
```

### **2. Remove indels/multi-allelic SNPs and apply missingness filters:**
```bash
plink --bfile $INPUT --snps-only just-acgt --make-bed --out ${INPUT}_snps_only
plink --bfile ${INPUT}_snps_only --geno 0.05 --make-bed --out ${INPUT}_geno_filtered
plink --bfile ${INPUT}_geno_filtered --mind 0.02 --make-bed --out ${INPUT}_mind_filtered
```
### **3.  Heterozygosity and missingness outlier detection:**
### Use Script heterozygosity_outliers.R in R

### **4.  Sex discrepancy check:**
```bash
plink --bfile $INPUT --split-x b37 --make-bed --out ${INPUT}_splitx
plink --bfile ${INPUT}_splitx --check-sex --out ${INPUT}_sexcheck
```
### Run R script to process sex check results
### Use Script plot_sexcheck.R in R

### **5.  Final QC filters (MAF, HWE):**
```bash
plink --bfile $INPUT --geno 0.02 --make-bed --out ${INPUT}_geno_final
plink --bfile ${INPUT}_geno_final --maf 0.01 --make-bed --out ${INPUT}_maf_filtered
plink --bfile ${INPUT}_maf_filtered --hwe 1e-10 --make-bed --out ${INPUT}_hwe_filtered
```

### **6.   Post-QC IBD/PI_HAT check:**
```bash
plink --bfile $INPUT --genome --min 0.2 --out ${INPUT}_final_ibd
sort -rk 10 ${INPUT}_final_ibd.genome > ${INPUT}_final_ibd_sorted.genome
```

---

# 6. LiftOver (only necessary if data is in Hg19)

![](https://github.com/ccmaues/pgt_images_github/blob/main/02_lift_over.png?raw=true)

A liftover is essential in genomics to convert genomic coordinates (e.g., gene positions, variants, or annotations) between different reference genome assemblies. Although the TOPMed imputation server accepts data in Hg19, this conversion step is crucial for Imputation QC, as the reference panel used in these steps is based on Hg38.

This step was performed with `LiftOverPlink` and `LiftOver`. Previous work ([Ormond et al., 2021](https://doi.org/10.1093/bib/bbab069) and [Github](https://github.com/cathaloruaidh/genomeBuildConversion)) has highlighted unusual behaviour in build conversion, such as SNVs mapping to a different chromosome names as conversion-unstable positions (CUPs), so these regions must be removed before LiftOver.

*Requirements to run:*

- [CUP files hg19](https://raw.githubusercontent.com/cathaloruaidh/genomeBuildConversion/master/CUP_FILES/FASTA_BED.ALL_GRCh37.novel_CUPs.bed)
- [LiftoverPlink](https://github.com/sritchie73/liftOverPlink)
- [Chain from hg19 to hg38](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)

*Main Command:*
```{bash}
awk '{print substr($1, 4), $2, $3, "SETID"}' "$liftover"/FASTA_BED.ALL_GRCh37.novel_CUPs.bed > "$liftover"/ALL_GRCh37_novel_CUPs_plink.bed

plink --bfile "$plink_file" --exclude range "$liftover"/ALL_GRCh37_novel_CUPs_plink.bed --output-chr M --make-bed --out "$plink_file"_del1
plink --bfile "$plink_file"_del1 --recode --output-chr M --out "$plink_file"
python3 "$liftover"/liftOverPlink.py --bin "$liftover"/liftOver --map "$plink_file".map --out "$path_file"/lifted --chain "$liftover"/hg19ToHg38.over.chain.gz
python3 "$liftover"/rmBadLifts.py --map "$path_file"/lifted.map --out "$path_file"/good_lifted.map --log "$path_file"/bad_lifted.dat
cut -f 2 "$path_file"/bad_lifted.dat > "$path_file"/to_exclude.dat
cut -f 4 "$path_file"/lifted.bed.unlifted | sed "/^#/d" >> "$path_file"/to_exclude.dat 
plink --file "$plink_file" --recode --out "$path_file"/lifted --exclude "$path_file"/to_exclude.dat 
plink --ped "$path_file"/lifted.ped --map "$path_file"/good_lifted.map --recode --out "$output"
plink --ped "$output".ped --map "$output".map --make-bed --out "$output"
```

---

# 7. Imputation QC:

Prior to imputation, we applied scripts and steps from [TOPMed's Imputation Server](https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/) to prepare the merged data for imputation. This included checking the data for compatibility with the TOPMed reference panel, filtering variants, and ensuring that the data was in the correct format for imputation.

<!---
Put parameters of imputation here
+ The liftover was done before the VCF creation, right?
--->

| Imputation service | Parameter |
| --- | --- |
| **Panel** | TOPMed r3 |
| **Array build** | GRCh38/hg38 |
| **rsq filter** | off |
| **Phasing** | Eagle v2.4 |
| **Population** | TOPMed Panel |
| **Mode** | Quality Control & Imputation |

*Requirements to run:*

- [HRC or 1000G Imputation preparation and checking](https://www.chg.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip)
- [TOPMed Reference Panel](https://legacy.bravo.sph.umich.edu/freeze5/hg38/download)
- [Script to format panel](https://www.chg.ox.ac.uk/~wrayner/tools/CreateTOPMed.zip)
- GRCh38 human genome reference

<!---
i lost genome link...?
--->

*Main Commands:*

```{bash}
plink --freq --bfile ${data} --out ${data}_freq
perl ${check}/HRC-1000G-check-bim.pl -b ${data}.bim -f ${data}_freq.frq -r ${check}/PASS.Variants.TOPMed_freeze5_hg38_dbSNP.tab.gz -h
bash Run-plink.sh
```

```{bash}
find ${data_path} -type f -name "*chr23*" | while IFS= read -r file; do
	new_name="${file//chr23/chrX}";
	mv "$file" "$new_name";
done
```

```{bash}
for i in {1..22} X; do
	plink --vcf ${data}-updated-chr${i}.vcf --output-chr chrMT --recode vcf --out ${data}-updated-chr${i} --noweb;
done
```

```{bash}
for i in {1..22} X; do
	bcftools +fixref ${data}-updated-chr${i}.vcf -Oz -o ${data}-updated-chr${i}.ref.vcf -- -d -f ${check}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -m flip 2>&1 | tee ${data}-updated-chr${i}.ref.vcf.log;
done
```

```{bash}
for i in {1..22} X; do vcf-sort ${data}-updated-chr${i}.ref.vcf | bgzip -c > ${output_path}/${batch}_chr${i}.vcf.gz; done
mv ${data}-updated* ${data_path}/*-HRC.txt Run-plink.sh ${data_prep}
```

*Removal*

| Dataset | liftover | Excluded | Unsolved | Allele mismatch | TOPMed QC | Final |
| -------- | --- | --- | --- | --- | --- | --- |
| **GSA** | 253 | 16,593 | 1,720 | 13 | 6,618 | X |
| **PSYCH** | 768 | 9,428 | 4,043 | 28 | 4,061 | X |
| **OMNI** | 824 | 19,046 | 1,437 | 48 | 10,800 | 639,434 |

<!---
put info later - might transfer to QC topic instead of imputation one
--->

*SNVs and samples per dataset*

| Dataset | GSA | PSYCH | OMNI |
| -------- | --- | --- | --- |
| **Initial Variants** | XXXX | XXXX | XXXX |
| **Final Variants** | 325,155 | 291,482 | 639,434 |
| | | | |
| **Initial samples** | XXX | XXX | XXX |
| **Final samples** | XXX | XXX | XXX |

---

# 8. Post imputation processing

![](https://raw.githubusercontent.com/ccmaues/pgt_images_github/refs/heads/main/08_post_impt.png)

The Rsq is a measure of the quality of imputed genotypes - higher values indicate a better quality - and the value is provided by `Minimac4` in the `TOPMed's imputation pipeline`. We set a threshold of 0.8 for the Rsq score (R²), which is a common threshold used in genetic studies to ensure high-quality imputed genotypes (e.g. **2. INFO filter** below).

*Main Commands:*

### **1. Unzip:**

```{bash}
for i in {1..22} X; do 
  7z e chr_"$i".zip -p'password' -o"$output"
done
```

### **2. INFO filter:**

```{bash}
find "$output_path" -name "*.dose.vcf.gz" | sort -V  | parallel -j 5 'bcftools view -i "R2>.8" -Oz -o {.}.filtered_r2_08.vcf.gz {}'
find "$output_path" -name "*.filtered_r2_08.vcf.gz" | sort -V | xargs bcftools concat -Oz -o $output_INFO
```

### **3. VCF to PLINK format:**

```{bash}
bcftools index -t "$output_INFO".vcf.gz
plink --vcf "$output_INFO".vcf.gz --make-bed --out "$output_INFO_PLINK"
```

<!---
put info later
--->

| Filter | INFO |
| --- | --- |
| **Psych** | X |
| **Omni** | X |
| **GSA** | X |

---

# 9. Quality control 3

Summary

<!---
put info later
--->

*Main Command:*

```{bash}
plink --bfile "$output_INFO_PLINK" --maf 0.01 --geno 0.02 -hwe 1e-10 --make-bed --out "$output_INFO_PLINK"_QC
plink --bfile "$output_INFO_PLINK"_QC --snps-only --make-bed --out "$output_INFO_PLINK"_QC_snp
plink2 --bfile "$output_INFO_PLINK"_merge_nodup_common --rm-dup exclude-all --make-bed --out "$output_INFO_PLINK"_merge_nodup_common_snp
```

*SNV removal at each step*

| Filter | Psych | Omni | GSA |
| --- | --- | --- | --- |
| **HWE** | 203,354 | 277,098 | 583,915 |
| **MAF** | 14,768,441 | 18,981,364 | 36,544,067 |
| **SNPs-Only** | 565247 | 588953 | 579855 |

*Final SNV number after QC*

| Filter | SNVs | Samples |
| --- | --- | --- |
| **Psych** | 9,884,542 | X |
| **Omni** | 10,324,705 | X |
| **GSA** | 10,138,583 | X |
---

# 10. Chip dataset merge

By merging datasets, we ensure that the final separated groups are derived from a unified, high-quality dataset, enabling robust downstream analyses. In this step, the three different chipsets were merged.

*Main Command:*

Refer to `4. Dataset merge` command description.

*Merged dataset*

| Dataset | SNVs | Samples |
| --- | --- | --- |
| Psych + Omni | 9,151,967 | 853 |
| Psych + Omni + Psych | 8,787,496 | 6,523 |

---

# 11. Cohort separation

Summary

*Install to run*

*External files to use:*

*Main Command:*

```{bash}

```

```{bash}

```

```{bash}

```

```
```
