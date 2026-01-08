# Array QC
# set parameters
inPath='/disk2/cgomes/08-2024-imputation-pipeline/roc_input/ROC_GDA-8v1-0_D1_hg37_multisample_gtc2vcf'
cleanedPath='/disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned'
qcstatsPath='/disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_stats'
clvcfPath='/disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned_vcfs'
clvcfbychrPath='/disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned_vcfs_bychrs'
refFile='/disk2/cgomes/08-2024-imputation-pipeline/roc_input/liftover/hg19_no_chr.fa.gz'

# cleaned plink data
# This steps will usew PLINK 1.9 and KING for QC and relatedness inference
plink --bfile ${inPath} --maf 0.01 --hwe 1e-6 --geno 0.05 --mind 0.05 --make-bed --out ${cleanedPath}

# estimate relationship
#king -b ${cleanedPath}.bed --degree 2 --related --prefix ${qcstatsPath}/MxGDAR_Nofilters_nodup_cleaned_king

# estimate sex
plink --bfile ${cleanedPath} --impute-sex --make-bed --out ${qcstatsPath}/MxGDAR_Nofilters_nodup_cleaned_sexplink

# estimate heterozygosity
plink --bfile ${cleanedPath} --het --out ${qcstatsPath}/MxGDAR_Nofilters_nodup_cleaned_hetplink


# Prepare data for imputation
# Fix headers - specific for ROC cohort
awk '{print $1 "-" $2 " " $2 " " $3 " " $4 " " $5 " " $6}' /disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned/plink_cleaned.fam > /disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned/temp_plink_cleaned.fam
rm /disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned/plink_cleaned.fam
mv /disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned/temp_plink_cleaned.fam /disk2/cgomes/08-2024-imputation-pipeline/roc_input/plink_cleaned/plink_cleaned.fam

# PLINK create vcf command
plink --bfile ${cleanedPath} --recode vcf-fid bgz --out ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned

# index with BCFtools
bcftools index ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned.vcf.gz

# convert vcf file to bcf filels
bcftools convert ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned.vcf.gz -Ob -o ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned.bcf.gz

# BCFtools run fix statistics       
export PATH="/usr/bin/bcftools-1.9:$PATH"
export BCFTOOLS_PLUGINS="/usr/bin/bcftools-1.9/plugins/"

bcftools +fixref ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned.bcf.gz -- -f ${refFile}

# BCFtools fix vcf reference command
bcftools +fixref ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned.bcf.gz -Oz -o ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned_fixed.vcf.gz -- -d -m flip -f ${refFile}

# Remove duplicared sites
bcftools norm -d all -f ${refFile} -Oz -o ${clvcfPath}/roc_forimputation_nodup.vcf.gz ${clvcfPath}/MxGDAR_Nofilters_nodup_cleaned_fixed.vcf.gz

# The above command might have changed the coordinates, we must sort the VCF.
bcftools sort ${clvcfPath}/roc_forimputation_nodup.vcf.gz -Oz -o ${clvcfPath}/roc_forimputation_sorted.vcf.gz

# Check VCF
python2 /disk2/cgomes/vcfcheck/checkVCF.py -r ${refFile} \
  -o ${clvcfPath}/roc_forimputation_sorted ${clvcfPath}/roc_forimputation_sorted.vcf.gz

# Iterate over chromosome numbers to transform vcf to compress vcf files
for chr_num in {1..22}
do
plink2 --vcf ${clvcfPath}/roc_forimputation_sorted.vcf.gz --chr ${chr_num} --recode vcf-4.2 bgz --out ${clvcfbychrPath}/roc_forimputation_chr${chr_num}
done

# Copy to s3
aws s3 sync --exclude="*" --include="*.vcf.gz" ${clvcfbychrPath} s3://prs-roc/raw/vcf/for_imputation/
