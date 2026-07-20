#!/bin/bash

dragen="<PATH_TO_DRAGEN_EXECUTABLE>"
cred_cfg="<PATH_TO_DRAGEN_LICENSE_CREDENTIALS>"
ref_fasta="<https://drive.google.com/drive/folders/11TCIk4XEgbzg5j9GRe_vvzaDEc1zFD11?usp=drive_link>"

config_dir="config"
output_dir="output"

number_of_batches=4
num_shards=102
num_threads=72

version=1
merge=1

###############################################################################
# Step 1: Aggregate gVCF files
#
# Loop:
#   batch = 1 to 10 (preferably around 1000 samples per batch as recommended in Illumina-style workflows)
#   shard = 1 to 102 (genome coordinate shard, where the whole genome is split into 102 shards)
#
# Input:
#   config/batch-gvcf/batch-1.gvcfs.list
#   config/batch-gvcf/batch-2.gvcfs.list
#   config/batch-gvcf/batch-3.gvcfs.list
#   .
#   .
#   .
#   config/batch-gvcf/batch-10.gvcfs.list
#
# Output:
#   output/step1-aggregate-gvcf/batch-*/shard-*/
###############################################################################

for batch in $(seq 1 10); do
    for shard in $(seq 1 102); do

        echo ""
        echo "STEP 1: Aggregate gVCF files"
        echo "Running batch ${batch}/10 and shard ${shard}/102"
        echo "Batch means approximately 1000 sample gVCFs."
        echo "Shard means one of 102 genome-coordinate splits."

        step_output_dir="${output_dir}/step1-aggregate-gvcf/batch-${batch}/shard-${shard}"

        "${dragen}" \
            --sw-mode \
            --lic-credentials "${cred_cfg}" \
            --enable-gvcf-genotyper-iterative true \
            --gvcfs-to-cohort-census true \
            --variant-list "${config_dir}/batch-gvcf/batch-${batch}.gvcfs.list" \
            --shard "${shard}/102" \
            --num-threads 72 \
            --logging-to-output-dir true \
            --output-directory "${step_output_dir}" \
            --ht-reference "${ref_fasta}" \
            --gg-min-qual 3.0 \
            --gg-remove-nonref true \
            --gg-discard-ac-zero true

    done
done

###############################################################################
# Step 2: Aggregate census files
#
# Loop:
#   shard = 1 to 102
#
# Input:
#   Step 1 census files:
#   output/step1-aggregate-gvcf/batch-*/shard-*/dragen.cns.gz
#
# Output:
#   output/step2-aggregate-census/version-1/shard-*/
###############################################################################

for shard in $(seq 1 102); do

    echo ""
    echo "STEP 2: Aggregate census files"
    echo "Running shard ${shard}/102"
    echo "This step aggregates census files across all batches for each shard."

    step_output_dir="${output_dir}/step2-aggregate-census/version-1/shard-${shard}"

    > "${step_output_dir}/census.list"

    for batch in $(seq 1 4); do
        realpath "${output_dir}/step1-aggregate-gvcf/batch-${batch}/shard-${shard}/dragen.cns.gz" >> "${step_output_dir}/census.list"
    done

    "${dragen}" \
        --sw-mode \
        --lic-credentials "${cred_cfg}" \
        --enable-gvcf-genotyper-iterative true \
        --aggregate-censuses true \
        --variant-list "${step_output_dir}/census.list" \
        --shard "${shard}/102" \
        --num-threads 72 \
        --logging-to-output-dir true \
        --output-directory "${step_output_dir}" \
        --ht-reference "${ref_fasta}" \
        --gg-max-subregion-size 1000 \
        --gg-concurrency-region-size 100000 \
        --verbose

done

###############################################################################
# Step 3: Generate multisample VCF
#
# Loop:
#   batch = 1 to 4
#   shard = 1 to 102
#
# Input:
#   Step 1 cohort file: dragen.cht.gz
#   Step 1 census file: dragen.cns.gz
#   Step 2 global census file: dragen.cns.gz
#
# Output:
#   output/step3-generate-msvcf/version-1/batch-*/shard-*/
###############################################################################

for batch in $(seq 1 4); do
    for shard in $(seq 1 102); do

        echo ""
        echo "STEP 3: Generate multisample VCF"
        echo "Running batch ${batch}/4 and shard ${shard}/102"
        echo "This step generates batch-level multisample VCFs for each shard."

        step_output_dir="${output_dir}/step3-generate-msvcf/version-1/batch-${batch}/shard-${shard}"

        "${dragen}" \
            --sw-mode \
            --lic-credentials "${cred_cfg}" \
            --enable-gvcf-genotyper-iterative true \
            --generate-msvcf true \
            --input-cohort-file "${output_dir}/step1-aggregate-gvcf/batch-${batch}/shard-${shard}/dragen.cht.gz" \
            --input-census-file "${output_dir}/step1-aggregate-gvcf/batch-${batch}/shard-${shard}/dragen.cns.gz" \
            --input-global-census-file "${output_dir}/step2-aggregate-census/version-1/shard-${shard}/dragen.cns.gz" \
            --shard "${shard}/102" \
            --num-threads 72 \
            --logging-to-output-dir true \
            --output-directory "${step_output_dir}" \
            --ht-reference "${ref_fasta}" \
            --gg-write-sample-qual true \
            --gg-max-subregion-size 1000 \
            --gg-concurrency-region-size 100000 \
            --verbose

    done
done

###############################################################################
# Step 4: Merge multisample VCFs per shard
#
# Loop:
#   shard = 1 to 102
#
# Input:
#   Step 3 VCF files from all batches for each shard
#
# Output:
#   output/step4-merge-shard-msvcf/version-1/merge-1/shard-*/
###############################################################################

for shard in $(seq 1 102); do

    echo ""
    echo "STEP 4: Merge multisample VCFs per shard"
    echo "Running shard ${shard}/102"
    echo "This step merges batch-level multisample VCFs for each shard."

    step_output_dir="${output_dir}/step4-merge-shard-msvcf/version-1/merge-1/shard-${shard}"


    > "${step_output_dir}/msvcf.list"

    for batch in $(seq 1 10); do
        realpath "${output_dir}/step3-generate-msvcf/version-1/batch-${batch}/shard-${shard}/dragen.vcf.gz" >> "${step_output_dir}/msvcf.list"
    done

    "${dragen}" \
        --sw-mode \
        --lic-credentials "${cred_cfg}" \
        --enable-gvcf-genotyper-iterative true \
        --merge-batches true \
        --gg-enable-indexing true \
        --input-batch-list "${step_output_dir}/msvcf.list" \
        --num-threads 72 \
        --logging-to-output-dir true \
        --output-directory "${step_output_dir}"

done

###############################################################################
# Step 5: Concatenate shard-level VCFs into chromosome-level VCFs per batch
#
# Loop:
#   batch = 1 to 10
#   chrom = 1 to 25
#
# Chromosome index:
#   1-22 = chr1-chr22
#   23   = chrX
#   24   = chrY
#   25   = chrM
#
# Input:
#   Step 3 shard-level VCF files
#
# Output:
#   output/step5-concat-batch-msvcf/version-1/batch-*/chrom-*/
#
# Note:
#   This step uses bcftools, not DRAGEN.
###############################################################################

module load bcftools-1.18

for batch in $(seq 1 10); do
    for chrom in $(seq 1 25); do

        echo ""
        echo "STEP 5: Concatenate shard-level VCFs into chromosome-level VCF"
        echo "Running batch ${batch}/4 and chromosome index ${chrom}/25"
        echo "Chromosome index 1-25 represents chr1-chr22, chrX, chrY, and chrM."

        step_output_dir="${output_dir}/step5-concat-batch-msvcf/version-1/batch-${batch}/chrom-${chrom}"

        rm -rf "${step_output_dir}"
        mkdir -p "${step_output_dir}"

        chrom_name=$(awk -v chrom="${chrom}" '$1==chrom {print $2}' "${config_dir}/chrom.25.list")

        > "${step_output_dir}/msvcf.list"

        awk -v chrom_name="${chrom_name}" '$2 ~ "^"chrom_name":" {print $1}' "${config_dir}/shards.hg38.3366.102.1.1.txt" | while read shard; do
            realpath "${output_dir}/step3-generate-msvcf/version-1/batch-${batch}/shard-${shard}/dragen.vcf.gz" >> "${step_output_dir}/msvcf.list"
        done

        bcftools concat \
            --naive-force \
            --no-version \
            -Oz \
            -o "${step_output_dir}/dragen.vcf.gz" \
            -f "${step_output_dir}/msvcf.list"

        bcftools index -t -f "${step_output_dir}/dragen.vcf.gz"

    done
done

###############################################################################
# Step 6: Merge chromosome-level multisample VCFs across batches
#
# Loop:
#   chrom = 1 to 25
#
# Input:
#   Step 5 chromosome-level VCF files from all batches
#
# Output:
#   output/step6-merge-chrom-msvcf/version-1/merge-1/chrom-*/
###############################################################################

for chrom in $(seq 1 25); do

    echo ""
    echo "STEP 6: Merge chromosome-level multisample VCFs across batches"
    echo "Running chromosome index ${chrom}/25"
    echo "This step creates the final merged chromosome-level VCF."

    step_output_dir="${output_dir}/step6-merge-chrom-msvcf/version-1/merge-1/chrom-${chrom}"

    rm -rf "${step_output_dir}"
    mkdir -p "${step_output_dir}"

    > "${step_output_dir}/msvcf.list"

    for batch in $(seq 1 10); do
        realpath "${output_dir}/step5-concat-batch-msvcf/version-1/batch-${batch}/chrom-${chrom}/dragen.vcf.gz" >> "${step_output_dir}/msvcf.list"
    done

    "${dragen}" \
        --sw-mode \
        --lic-credentials "${cred_cfg}" \
        --enable-gvcf-genotyper-iterative true \
        --merge-batches true \
        --gg-enable-indexing true \
        --input-batch-list "${step_output_dir}/msvcf.list" \
        --num-threads 72 \
        --logging-to-output-dir true \
        --output-directory "${step_output_dir}"

done

echo ""
echo "All 6 IGG joint-calling steps completed."