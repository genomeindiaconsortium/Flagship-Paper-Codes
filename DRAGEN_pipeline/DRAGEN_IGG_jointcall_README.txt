# DRAGEN / IGG Joint Calling Workflow

Illumina dragen documentation: https://support-docs.illumina.com/SW/DRAGEN_v310/Content/SW/DRAGEN/gVCFGenotyper.htm
The workflow starts from individual sample gVCF files and produces cohort-level multisample VCF files.

## Important note

This repository does not include Illumina DRAGEN software license files or credentials. Users must obtain the required DRAGEN software and license directly from Illumina.

## Required inputs

Before running the script, prepare:

- Per-sample gVCF files
- One gVCF list file per batch
- Reference FASTA file < https://drive.google.com/drive/folders/11TCIk4XEgbzg5j9GRe_vvzaDEc1zFD11?usp=drive_link>
- DRAGEN license credentials
- DRAGEN executable
- Chromosome mapping file
- Shard mapping file

Example gVCF list file:

```text
/path/to/sample1.g.vcf.gz
/path/to/sample2.g.vcf.gz
/path/to/sample3.g.vcf.gz

batch = group of sample gVCFs, preferably around 1000 samples per batch
shard = genome-coordinate split; the genome is divided into 102 shards
chrom = chromosome index; 1-22 = chr1-chr22, 23 = chrX, 24 = chrY, 25 = chrM

Step 1: Aggregate gVCF files
        Runs for each batch and each shard.

Step 2: Aggregate census files
        Runs for each shard.

Step 3: Generate multisample VCF files
        Runs for each batch and each shard.

Step 4: Merge multisample VCFs per shard
        Runs for each shard.

Step 5: Concatenate shard-level VCFs into chromosome-level VCFs
        Runs for each batch and each chromosome.
        This step uses bcftools.

Step 6: Merge chromosome-level VCFs across batches
        Runs for each chromosome.