#!/bin/bash
##### Script to compute PGS for externalizing phenotype using PRSCs #####
#### Catarina, July 2025 ####

set -e

# Step 0 - Perform liftover if necessary
# Step 1 - Estimate posterior effect sizes using PRS-CS
# Step 2 - Compute individual level PGS using PLINK

######################## Variables ####################
wd="/disk2/cgomes/2025-12-04-pgs-ext"
ORIGINAL_SUMSTATS="$wd/sumstas/GSEM.GWAS.EXTERNALIZING.excl23andMe.adjusted.txt"
GWAS_SAMPLE_SIZE=1045957
POPULATION="EUR"
GWAS_BUILD="hg19"  # Set to "hg38" to skip liftover - roc is in hg38

# LD matrix is in GRCh37
PATH_TO_REFERENCE="/disk2/datasets/ldmatrix/ldblk_1kg_eur" # Downloaded from https://github.com/getian107/PRScs
OUTPUT_DIR="/disk2/cgomes/2025-12-04-pgs-ext/output"
OUTPUT_FILE_PREFIX="gwas_externalizing"

VCF_INPUT="/disk2/cgomes/2024-23-12-post-imputation-qc/ROC_GDA.imputed.fixed.filtered.dose.vcf.gz"

VCF_LIFTED="$wd/lifted/ROC_GDA.imputed.fixed.filtered.dose.hg19.vcf.gz"
PLINK_PREFIX_HG19="$wd/plink/ROC_GDA.imputed.fixed.filtered.dose.hg19"
PLINK_PREFIX_HG38="$wd/plink/ROC_GDA.imputed.fixed.filtered.dose.hg38"

# Detect system resources
NUM_THREADS=$(nproc)
JAVA_MEM=$(awk '/MemTotal/ {printf "%.0f", $2/1024/1024 * 0.8}' /proc/meminfo)  # 80% of total RAM in GB

# Make output dirs
#mkdir -p $OUTPUT_DIR/prs $wd/lifted $wd/plink

#########################################
### Step 0 - Liftover VCF if needed
#########################################

if [[ "$GWAS_BUILD" == "hg19" ]]; then

    echo "Renaming VCF contigs to match reference for liftover..."
    
    bcftools annotate \ --threads $NUM_THREADS \
    --rename-chrs <(for i in {1..22} X Y MT; do echo -e "$i\tchr$i"; done) \
    $VCF_INPUT -Oz -o $wd/lifted/ROC_GDA.with_chr.vcf.gz

    tabix --threads $NUM_THREADS -p vcf $wd/lifted/ROC_GDA.with_chr.vcf.gz
    VCF_INPUT="$wd/lifted/ROC_GDA.with_chr.vcf.gz"

    echo "GWAS build is hg19. Lifting over VCF from hg38 to hg19..."

    # Run Picard LiftoverVcf
    
    picard -Xmx${JAVA_MEM}g LiftoverVcf \
        I=$VCF_INPUT \
        O=$VCF_LIFTED \
        CHAIN=/disk2/datasets/hg38ToHg19.over.chain \
        REJECT=$wd/lifted/rejected_variants.vcf.gz \
        R=/disk2/datasets/hg19.fa \
        MAX_RECORDS_IN_RAM=100000

    echo "Removing non-canonic chr"
    CANON_VCF="$wd/lifted/ROC_GDA.hg19.canonical.vcf.gz"

    bcftools view --threads "$NUM_THREADS" \
    -r 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM' \
    -Oz -o "$CANON_VCF" "$VCF_LIFTED"
    
    tabix --threads "$NUM_THREADS" -p vcf "$CANON_VCF"

    echo "Converting lifted VCF to PLINK format..."
    plink \
      --vcf $CANON_VCF \
      --chr 1-22 \
      --double-id \
      --out $PLINK_PREFIX_HG19 \
      --memory 20000 \
      --threads $NUM_THREADS

    VALIDATION_BIM_PREFIX=$PLINK_PREFIX_HG19

    echo "Removing chr prefix"
    sed -i 's/^chr//g' "$PLINK_PREFIX_HG19.bim"

    echo "create a de-duplicated bed/bim/fam..."
    
    # create a de-duplicated bed/bim/fam
    plink2 \
    --bfile "${VALIDATION_BIM_PREFIX}" \
    --set-missing-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 200 truncate \
    --make-bed \
    --out "${VALIDATION_BIM_PREFIX}_idfix" \
    --threads "${NUM_THREADS}" \
    --memory 20000

    plink2 \
    --bfile "${VALIDATION_BIM_PREFIX}_idfix" \
    --rm-dup exclude-all \
    --make-bed \
    --out "${VALIDATION_BIM_PREFIX}_clean" \
    --threads "${NUM_THREADS}" \
    --memory 20000
    
else
    echo "GWAS build is hg38. No liftover needed — converting original VCF..."
    plink \
      --vcf $VCF_INPUT \
      --id-delim \
      --out $PLINK_PREFIX_HG38 \
      --memory 20000

    VALIDATION_BIM_PREFIX=$PLINK_PREFIX_HG38
fi


#########################################
### Step 1 - Run PRS-CS
#########################################

echo "Running PRS-CS..."

python /home/ubuntu/PRScs/PRScs.py \
 --ref_dir=$PATH_TO_REFERENCE \
 --bim_prefix=$VALIDATION_BIM_PREFIX \
 --n_gwas=$GWAS_SAMPLE_SIZE \
 --sst_file=$ORIGINAL_SUMSTATS \
 --out_dir=$OUTPUT_DIR \
 --out_name=$OUTPUT_FILE_PREFIX

### Combine per-chr posterior files ###
POSTERIOR="${OUTPUT_DIR}/${OUTPUT_FILE_PREFIX}_pst_eff_a1_b0.5_phiauto_all_chr.txt"
> "$POSTERIOR"

for chr_num in {1..22}; do
    file="${OUTPUT_DIR}/output_pst_eff_a1_b0.5_phiauto_chr${chr_num}.txt"    
    cat "$file" >> "$POSTERIOR"
done

#########################################
### Step 2 - Compute individual PGS
#########################################

echo "Computing individual-level PGS..."
plink \
--score $POSTERIOR 2 4 6 sum center \
--bfile ${VALIDATION_BIM_PREFIX}_clean \
--out ${OUTPUT_DIR}/prs/ROC_GDA_externalizing.${GWAS_BUILD} \
--allow-no-sex \
--memory 20000

echo "Pipeline complete ✅"
