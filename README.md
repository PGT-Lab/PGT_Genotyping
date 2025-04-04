# PGT quality control pipeline

### Processing of the genotping array datasets at Laboratory of Integrative Neurosciences (LINC) at UNIFESP.


<!---
put i here a figure showing the "steps" done
<img src="https://cdn-icons-png.flaticon.com/512/8704/8704721.png" width="35" height="35">
--->

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
| Lucas Toshio Ito | XX | XX |
| Adrielle | XX | XX |
| Gustavo | XX | XX |

---

# 1. IDAT to PLINK

![](https://github.com/ccmaues/pgt_images_github/blob/main/01_idat_to_vcf.png?raw=true)

This script convertes IDAT files to GTC format, then to VCF format, and finally to PLINK format for downstream analysis. It uses `bcftools` plugins for file conversions, reheaders and normalizes the VCF file, indexes it for efficient querying, and converts it to PLINK format.

*Requirements for this step:*

- [gtc2vcf](https://github.com/freeseek/gtc2vcf)
- [Chip manifest and cluster files](https://support.illumina.com/array/downloads.html)
- [GRCh37 human genome reference](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz)
- [GRCh38 human genome reference](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)

```{bash}
sudo apt install wget unzip git g++ zlib1g-dev bwa unzip samtools msitools cabextract mono-devel libgdiplus icu-devtools bcftools
sudo apt install libbz2-dev libssl-dev liblzma-dev libgsl0-dev
```

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
  bcftools norm \
  --threads 15 \
  --no-version \
  -Oz \
  -c x \
  -f "$ref" | \
  tee "$out_prefix.vcf.gz" | \
  bcftools index \
  --threads 15 \
  -ft \
  --output  "$out_prefix.vcf.gz.tbi"
```

### **3. VCF to PLINK v1.9 files (bed bim fam)**

```{bash}
plink --vcf "$out_prefix" --const-fid --update-sex "$sex_file" --make-bed --out "$out_prefix"
```

---

# 2. LiftOver

![](https://github.com/ccmaues/pgt_images_github/blob/main/02_lift_over.png?raw=true)
<!---
put i here a figure showing the "steps" done
--->

A liftover is essential in genomics to convert genomic coordinates (e.g., gene positions, variants, or annotations) between different reference genome assemblies. This step was performed with `LiftOverPlink` and `LiftOver`.

*Requirements to run:*

- [LiftoverPlink](https://github.com/sritchie73/liftOverPlink.git)
- [Chain from hg19 to hg38](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz)

*Main Command:*
```{bash}
plink --bfile "$plink_file" --recode --out "$prefix_hg19"
python3 liftOverPlink.py --bin ./liftOver --map "$prefix_hg19.map" --out lifted --chain hg19ToHg38.over.chain.gz
python3 rmBadLifts.py --map lifted.map --out good_lifted.map --log bad_lifted.dat
cut -f 2 bad_lifted.dat > to_exclude.dat
cut -f 4 lifted.bed.unlifted | sed "/^#/d" >> to_exclude.dat 
plink --file "$prefix_hg19" --recode --out lifted --exclude to_exclude.dat 
plink --ped lifted.ped --map good_lifted.map --recode --out "$output"
plink --ped "$output.ped" --map "$output.map" --make-bed --out "$output"
```

---

# 3. Check overlap

We estimated the SNV overlap between different chips and cohorts providing a summary of overlaps across those datasets.

![SNVs per chip version](https://github.com/ccmaues/pgt_images_github/blob/main/chip_variants.png?raw=true)

The table below shows the number of overlapping SNVs between different datasets. The diagonal cells represent the dataset itself, while the other cells show the number of overlapping SNVs between the corresponding datasets.

| Dataset                  | PC_Kings_2016_2017 | PI_ESALQ_2020 | PSC_Kings_2019 | PSC_Kings_2016_2019 | TP_USP_2020 | SC_CHOP_2017 | PS_CHOP_2017 |
|--------------------------|--------------------|---------------|---------------|---------------------|-------------|-------------|-------------|
| **PC_Kings_2016_2017**   | -                  | 109542        | 108461        | 587111              | 107784      | 258365      | 256052      |
| **PI_ESALQ_2020**        | 109542             | -             | 635892        | 109542              | 613312      | 139140      | 138358      |
| **PSC_Kings_2019**       | 108461             | 635892        | -             | 108461              | 589402      | 137739      | 136991      |
| **PSC_Kings_2016_2019**  | 587111             | 109542        | 108461        | -                   | 107784      | 258365      | 256052      |
| **TP_USP_2020**         | 107784             | 613312        | 589402        | 107784              | -           | 133277      | 132530      |
| **SC_CHOP_2017**         | 258365             | 139140        | 137739        | 258365              | 133277      | -           | 719023      |
| **PS_CHOP_2017**         | 256052             | 138358        | 136991        | 256052              | 132530      | 719023      | -           |

*Requirements to run:*

```{r}
install.package("data.table")
install.package("plyr")
```

---

# 4. Quality control 1

![](https://raw.githubusercontent.com/ccmaues/pgt_images_github/refs/heads/main/05_QC.png)

Quality control (QC) is critical in genetic studies to ensure data reliability by removing poor-quality samples and variants. This QC step was applied to avoid **false associations, population stratification bias, or spurious results**, compromising study validity. T the parameters values were the same as in RICOPILI pipeline, with the addition of relatedness filter.

<!---
add refs here for QC
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
| PC_Kings_2016_2017 | 257,032 | 7,612 | 3,609 | 0 | 0 | 0 | 0 | 0 | 0 |
| PI_ESALQ_2020 | 145,917 | 15,852 | 71,940 | 6 | 0 | 0 | 0 | 0 | 0 |
| PSC_Kings_2019 | 15,293 | 6,003 | 298,548 | 6 | 2,441| 0 | 0 | 0 | 0 |
| PSC_Kings_2016_2019 | 85,573 | 54,891 | 97,262 | 6 | 1,815 | 0 | 0 | 0 | 0 |
| TP_USP_2020 | 126,363 | 3 | 35,434 | 0 | 95 | 0 | 0 | 0 | 0 |
| SC_CHOP_2017 | 21,622 | 16,066 | 9,932 | 2 | 0 | 0 | 0 | 0 | 0 |
| PS_CHOP_2017 | 21,096 | 756 | 15,491 | 0 | 0 | 0 | 0 | 0 | 0 |

---

# 5. Dataset merge

The merged was performed only with the inteserction of SNVs in pairs of datasets.

*Main Command:*

### **1. SNV list:**

```{bash}
plink --bfile "$dataset1" --list-duplicate-vars --out "$dataset1_new"
plink --bfile "$dataset1" --exclude "$dataset1_new".dupvar --make-bed --out "$dataset1_new"_nodup
plink --bfile "$dataset1_new"_nodup --write-snplist --out "$dataset1_new"_nodup
sort -o "$dataset1_new"_nodup.snplist "$dataset1_new"_nodup.snplist
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
plink --bfile "$dataset1_new"_nodup_common --bmerge "$dataset2_new"_nodup_common_fid1 --make-bed --out "$psych_merged"
```

---

# 6. Imputation QC:

Prior to imputation, we applied scripts and steps from [TOPMed's Imputation Server](https://topmedimpute.readthedocs.io/en/latest/prepare-your-data/) to prepare the merged data for imputation. This included checking the data for compatibility with the HRC reference panel, filtering variants, and ensuring that the data was in the correct format for imputation.

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

*External files to use:*
```{bash}
wget https://www.chg.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.2.13-NoReadKey.zip
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
```

*Main Command:*

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
summary_file="unresolved_variants_summary.txt"
echo -e "Chromosome\tUnresolved_Variants" > ${data_path}/${summary_file}
total_unresolved=0
for i in {1..22} X; do
  log_file="${data}-updated-chr${i}.ref.vcf.log"
  if [ -f "$log_file" ]; then
    unresolved_count=$(grep -oP "^NS\s+unresolved\s+\K\d+" "$log_file")
    unresolved_count=${unresolved_count:-0}  # Default to 0 if empty
    total_unresolved=$((total_unresolved + unresolved_count))
    echo -e "chr${i}\t${unresolved_count}" >> ${data_path}/${summary_file}
  else
    echo -e "chr${i}\tN/A" >> ${data_path}/${summary_file}
  fi
done
echo -e "Total\t${total_unresolved}" >> ${data_path}/${summary_file}
cat ${data_path}/${summary_file}
# Sorting and zipping final file
for i in {1..22} X; do vcf-sort ${data}-updated-chr${i}.ref.vcf | bgzip -c > ${output_path}/${batch}_chr${i}.vcf.gz; done
# Moving all intermediary files to specific folder
mv ${data}-updated* ${data_path}/*-HRC.txt Run-plink.sh ${data_prep}
```

```{bash}
for i in {1..22} X; do
  7z e ${output_path_Omni}/chr_"$i".zip -p'OUad:y0CdWXx5W' -o${output_path_Omni}/unzip
done
```

*Removal*
| Dataset | liftover | Excluded | Unsolved | Allele mismatch | TOPMed QC | Final |
| -------- | --- | --- | --- | --- | --- | --- |
| GSA | 253 | 16593 | 1720 | 13 | 6618 |
| PSYCH | 768 | 9428 | 4043 | 28 | 4061 |
| OMNI | 824 | 19046 | 1437 | 48 | 10800 | 639434 |

<!---
put info later
--->

*SNVs and samples per dataset*

| Dataset | GSA | PSYCH | OMNI |
| -------- | --- | --- | --- |
| Initial Variants | XXXX | XXXX | XXXX |
| Final Variants | 325155 | 291482 | 639434 |
| | | | |
| Initial samples | XXX | XXX | XXX |
| Final samples | XXX | XXX | XXX |

---

# 7. Post imputation processing

The Rsq is a measure of the quality of imputed genotypes, with higher values indicating better quality. The value is provided by `Eagle` in the `TOPMed's imputation pipeline`. We set a threshold of 0.8 for the Rsq score, which is a common threshold used in genetic studies to ensure high-quality imputed genotypes.

<!---
maybe a ref here...? + check eagle output
--->

*Main Command:*

```{bash}
find "${output_path_GSA}" -name "*.dose.vcf.gz" |\
sort -V |\
parallel -j 5 'bcftools view -i "R2>.8" -Oz -o {.}.filtered_r2_08.vcf.gz {}'

find "${output_path_GSA}" -name "*.filtered_r2_08.vcf.gz" | sort -V | xargs bcftools concat -Oz -o "${output_path_GSA}/imputed_GSA_r2_08.vcf.gz"
```

```{bash}
bcftools index -t "${output_path_GSA}/imputed_GSA_r2_08.vcf.gz"
plink --vcf ${output_path_GSA}/imputed_GSA_r2_08.vcf.gz --make-bed --out ${output_path_GSA}/imputed_GSA_r2_08
```

---

# 8. Quality control 2

Summary

*Install to run*

*External files to use:*

*Main Command:*

```{bash}

```

---

# 9. Dataset merge

Summary

*Install to run*

*External files to use:*

*Main Command:*

```{bash}

```

---

# 10. Cohort separation

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