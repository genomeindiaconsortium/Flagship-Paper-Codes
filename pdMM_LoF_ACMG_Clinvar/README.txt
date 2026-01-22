README: 

pDMM_LoF_ACMG.py

Identification of putative Deleterious Missense (pDMM), High-Confidence Loss-of-Function (HC-LoF) Variants and filtering of ACMG Gene variants
This script processes variant annotation files and VCF data to identify two classes of potentially deleterious variants: putative deleterious missense mutations (pDMM) and high-confidence loss-of-function (HC-LoF) variants. It also cross-references against the ACMG secondary findings gene list to prioritize clinically actionable findings.

Key Analysis Steps:

pDMM Selection: Filters for MODERATE-impact (missense) variants predicted to be deleterious using REVEL ³ 0.644 and CADD_PHRED ³ 20 scores.
Homozygous HC-LoF Analysis:
Selects High-Confidence (HC) Loss-of-Function (LoF) variants(LOFTEE).
Extracts genotypes of these variants from the source VCF using bcftools.
Identifies variants with at least one homozygous alternate individual in the cohort.
ACMG Gene Filtering: Restricts both the pDMM and HC-LoF variant sets to genes listed in the ACMG recommendations (e.g., BRCA1, LDLR, TP53).

Outputs:

deleterious_variants.csv: A list of all missense variants meeting the deleterious scoring criteria (REVEL ³ 0.644 & CADD_PHRED ³ 20).
hc_lof_variants.csv: A list of all High-Confidence LoF variants identified in the annotation file.
hc_lof_homozygous_variants.csv: A subset of HC-LoF variants that are found in a homozygous state in at least one individual, including sample IDs.
acmg_pdmm_variants.csv: Deleterious missense variants restricted to ACMG genes.
acmg_hc_lof_variants.csv: High-Confidence LoF variants restricted to ACMG genes (regardless of zygosity).

Dependencies:

pandas (Python library)
bcftools (Command-line tool, must be in system PATH)
Input VCF file (input.vcf.gz)
Variant annotation file (CSV format)
ACMG gene list (acmg_gene_chrom.csv)


------------------------------------------------------------------------





Clinvar.py

Extraction of Population-Specific Minor Allele Frequencies (MAF) for >=2 star review status(>=2*) ClinVar Variants in GI dataset.

This utility script reads population frequency files of the 83 populations(.frq) of GI to extract ClinVar >=2* variants specific allele frequency data.




Key Analysis Steps:
ClinVar Filtering:
Loads the variant_summary_all.txt file.
Filters for GRCh38 assembly.
Retains only variants with >=2* review statuses: "criteria provided, multiple submitters, no conflicts", "reviewed by expert panel", or "practice guideline".
Population Frequency Matching:
Iterates through all 83 populations *.frq files in the directory (standard PLINK output format).
Parses genomic coordinates directly from SNP IDs.
Matches variants against the filtered ClinVar list (checking both Ref/Alt and Alt/Ref orientations).
Chunked Processing: Reads frequency files in chunks of 100,000 rows to optimize memory usage.

Outputs:
maf_intermediate/maf_{population}.tsv: Individual tab-separated files for each population containing:

Chromosome
Position
Reference Allele
Alternate Allele
Minor Allele Frequency (MAF)

Dependencies:
variant_summary_all.txt (ClinVar summary file downloaded on February 2025)
*.frq files (PLINK frequency files of 83 populations present in the working directory)

