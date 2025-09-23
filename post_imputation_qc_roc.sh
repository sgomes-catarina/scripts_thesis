################################################################################
################# ROC - QC of array data post imputation #######################
################################################################################

# This script runs QC on the ROC GDA data after TOPMed2 imputation 
# Script based on: https://github.com/jr-baldwin/ACEs_mental_health_RR/blob/main/3_ABCD_QC_genetic_20220204.sh

# ======= Load file paths ==============
# Load variables / directory paths
#childGeneticData="s3://prs-roc/raw/vcf/imputed/"

# Change directory to my file
cd /disk2/cgomes/2024-23-12-post-imputation-qc/

# =========== Download and unzip VCF files ================
aws s3 sync $childGeneticData . --exclude "*" --include "chr*"

# Unzip the info files
for i in {1..22}; do
gunzip chr_${i}/chr${i}.info.gz; 
done 

# =========== QC plots and metrics ================ #

Rscript post_imputation_qc_roc.R


# ======= Extract SNPs with imputation R2 > 0.8 with bcftools ======= #
# Cria um arquivo temporário com a lista de VCFs
vcf_list="vcfs_to_concat.txt"
> $vcf_list  # limpa o conteúdo caso já exista

for i in {1..22}; do
  echo "==> Processing chromosome $i..."

  input_vcf="chr_${i}/chr${i}.dose.vcf.gz"
  output_vcf="chr_${i}/chr${i}.filtered.dose.vcf.gz"

  bcftools view -i 'INFO/R2>=0.8' -Oz -o $output_vcf $input_vcf

  # Index the output VCF
  bcftools index $output_vcf

  echo "chr_${i}/chr${i}.filtered.dose.vcf.gz" >> $vcf_list
done

# Concatenate all filtered chromosomes
echo "==> Merging all filtered chromosomes into a single VCF..."
bcftools concat -Oz -o ROC_GDA.imputed.filtered.dose.vcf.gz -f $vcf_list

# Index the final VCF
bcftools index ROC_GDA.imputed.filtered.dose.vcf.gz

echo "✅ Done!"
