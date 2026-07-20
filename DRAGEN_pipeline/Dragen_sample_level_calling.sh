#!/bin/bash
set -euo pipefail


SAMPLE_NAME=$1
REF_DIR=$2
OUT_BASE_DIR=${3:-"/path/to/your/output"}
ADAPTER_R1=${4:-"illumina_adapter_r1.txt"}
ADAPTER_R2=${5:-"illumina_adapter_r2.txt"}

FASTQ_LIST="${SAMPLE_NAME}_FastqList.csv"
SAMPLE_OUT_DIR="${OUT_BASE_DIR}/${SAMPLE_NAME}"

mkdir -p "${SAMPLE_OUT_DIR}"

dragen -f \
    -r "${REF_DIR}" \
    --fastq-list "${FASTQ_LIST}" \
    --enable-variant-caller true \
    --vc-emit-ref-confidence GVCF \
    --output-directory "${SAMPLE_OUT_DIR}" \
    --trim-min-quality 15 \
    --trim-adapter-read1 "${ADAPTER_R1}" \
    --trim-adapter-read2 "${ADAPTER_R2}" \
    --output-file-prefix "${SAMPLE_NAME}" \
    --enable-duplicate-marking true \
    --enable-map-align-output true \
    --read-trimmers adapter,quality \
    --trim-min-length 50 \
    --enable-bam-indexing true

echo "DRAGEN pipeline completed for ${SAMPLE_NAME}."
