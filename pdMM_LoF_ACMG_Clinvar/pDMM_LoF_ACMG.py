
####--------------------------------pDMM criteria--------------------------------------
import pandas as pd

# ==== Step 1: Load annotation file ====
df = pd.read_csv("your_annotation_file.csv") 

# ==== Step 2: Select missense variants ====
# Only keep variants annotated as MODERATE impact (missense)
missense_variants = df[df['IMPACT'] == 'MODERATE'].copy()

# ==== Step 3: Apply deleterious criteria ====
# Ensure REVEL and PHRED are numeric
missense_variants['REVEL'] = pd.to_numeric(missense_variants['REVEL'], errors='coerce')
missense_variants['PHRED'] = pd.to_numeric(missense_variants['PHRED'], errors='coerce')

# Filter for deleterious missense variants
deleterious_variants = missense_variants[
    (missense_variants['REVEL'] >= 0.644) & 
    (missense_variants['PHRED'] >= 20)
]

# ==== Step 4: Save results ====
deleterious_variants.to_csv("deleterious_variants.csv", index=False)

print("Deleterious missense variants selection complete!")
print(f"Total deleterious variants: {deleterious_variants.shape[0]}")


###----------------------------------------homozygous HC-LoF criteria--------------------------------
Steps:
1. Select HC LoF variants from annotation file
2. Extract genotypes for these variants from a VCF
3. Count zygosity including homozygous ALT samples
4. Merge annotation with genotype counts and keep variants with ≥1 homozygous individual


import pandas as pd
import subprocess
import tempfile
import os

# -----------------------------
# Step 1: Select High-Confidence LoF variants
# -----------------------------
annot_file = "annotation_file.csv"  # Replace with your annotation file
df = pd.read_csv(annot_file, low_memory=False)

hc_lof = df[(df["LoF"] == "HC") & (df["PICK"] == 1)]
hc_lof.to_csv("hc_lof_variants.csv", index=False)
print(f"HC LoF variants selected: {hc_lof.shape[0]}")

# -----------------------------
# Step 2: Extract genotypes for HC LoF variants from VCF
# -----------------------------
vcf_file = "input.vcf.gz"  # Replace with your VCF file
variants = pd.read_csv("hc_lof_variants.csv", usecols=["Chr", "POS"])
variants["Chr"] = variants["Chr"].astype(str).str.replace("chr", "")
variants["POS"] = variants["POS"].astype(int)

with tempfile.NamedTemporaryFile("w", delete=False) as f:
    variants.to_csv(f.name, sep="\t", index=False, header=False)
    pos_file = f.name

# Get sample IDs from VCF
samples = subprocess.check_output(f"bcftools query -l {vcf_file}", shell=True, text=True).strip().split("\n")

out_file = "hc_lof_genotypes.txt"
with open(out_file, "w") as o:
    o.write("\t".join(["CHROM", "POS", "INFO"] + samples) + "\n")

# Extract variant genotypes
view_cmd = ["bcftools", "view", "-R", pos_file, vcf_file]
query_cmd = ["bcftools", "query", "-f", "%CHROM\t%POS\t%INFO[\t%GT]\n"]

with open(out_file, "a") as o:
    p1 = subprocess.Popen(view_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(query_cmd, stdin=p1.stdout, stdout=o)
    p1.stdout.close()
    p2.communicate()

os.remove(pos_file)
print("HC LoF genotypes extracted.")

# -----------------------------
# Step 3: Count homozygous LoF genotypes
# -----------------------------
def parse_genotype_counts(geno_file):
    genotype_counts = {}
    with open(geno_file, 'r') as f:
        header = f.readline().strip().split()
        sample_ids = header[3:]  # Samples start after CHROM, POS, INFO

        for line in f:
            parts = line.strip().split()
            key = f"{parts[0]}_{parts[1]}"
            genotypes = parts[3:]

            if key not in genotype_counts:
                genotype_counts[key] = {
                    'CHROM': parts[0],
                    'POS': parts[1],
                    'INFO': parts[2],
                    '0|0': 0,
                    '0|1': 0,
                    '1|0': 0,
                    '1|1': 0,
                    '.|.': 0,
                    'HOMO_SAMPLE_IDS': []
                }

            counts = genotype_counts[key]
            for i, gt in enumerate(genotypes):
                if gt in counts:
                    counts[gt] += 1
                    if gt == "1|1":
                        counts['HOMO_SAMPLE_IDS'].append(sample_ids[i])

    return genotype_counts

def genotype_counts_to_df(genotype_counts):
    rows = []
    for key, c in genotype_counts.items():
        rows.append({
            'CHROM': c['CHROM'],
            'POS': c['POS'],
            'INFO': c['INFO'],
            'HOM_REF': c['0|0'],
            'HET_01': c['0|1'],
            'HET_10': c['1|0'],
            'HOM_ALT': c['1|1'],
            'MISSING': c['.|.'],
            'HOMO_SAMPLE_IDS': ",".join(c['HOMO_SAMPLE_IDS'])
        })
    return pd.DataFrame(rows)

geno_counts = parse_genotype_counts(out_file)
geno_df = genotype_counts_to_df(geno_counts)
geno_df.to_csv("genotype_counts.txt", index=False, sep="\t")
print(f"Genotype counts computed for {geno_df.shape[0]} HC LoF variants.")

# -----------------------------
# Step 4: Merge annotation with genotype counts, keep variants with ≥1 homozygote
# -----------------------------
hc_lof["CHROM"] = hc_lof["Chr"].astype(str).str.replace("chr", "")

final_df = hc_lof.merge(geno_df, on=["CHROM", "POS"], how="inner")
final_df = final_df[final_df["HOM_ALT"] >= 1]
final_df.to_csv("hc_lof_homozygous_variants.csv", index=False)

print(f"HC LoF variants with ≥1 homozygous individual: {final_df.shape[0]}")


####--------------------------------ACMG Gene Filtering (pDMM & HC-LoF)--------------------------------
# This section is executed to filter variants restricted to genes found in the ACMG list.
# 1. pDMM variants (deleterious missense) are filtered.
# 2. HC LoF variants are filtered (regardless of zygosity).

# ==== Step A: The ACMG gene list is loaded ====
acmg_file = "acmg_gene_chrom.csv"
acmg_df = pd.read_csv(acmg_file)

# A set is created for efficient gene lookup
acmg_genes = set(acmg_df['Gene'])
print(f"ACMG gene list loaded. Total genes: {len(acmg_genes)}")

# ==== Step B: pDMM variants are filtered (if available) ====
if 'deleterious_variants' in locals():
    # The gene column is identified
    gene_col_missense = 'Gene'
    if 'Gene' not in deleterious_variants.columns and 'SYMBOL' in deleterious_variants.columns:
        gene_col_missense = 'SYMBOL'
    
    # Only deleterious missense variants within ACMG genes are retained
    acmg_pdmm = deleterious_variants[deleterious_variants[gene_col_missense].isin(acmg_genes)].copy()
    
    # The pDMM results are saved
    acmg_pdmm.to_csv("acmg_pdmm_variants.csv", index=False)
    print(f"pDMM filtering complete. Variants retained: {acmg_pdmm.shape[0]}")
else:
    print("pDMM filtering skipped: 'deleterious_variants' DataFrame was not found.")

# ==== Step C: HC LoF variants are filtered (if available) ====
if 'hc_lof' in locals():
    # The gene column is identified
    gene_col_lof = 'Gene'
    if 'Gene' not in hc_lof.columns and 'SYMBOL' in hc_lof.columns:
        gene_col_lof = 'SYMBOL'
    
    # Only HC LoF variants within ACMG genes are retained (Homozygosity check is bypassed)
    acmg_lof_variants = hc_lof[hc_lof[gene_col_lof].isin(acmg_genes)].copy()

    # The HC LoF results are saved
    acmg_lof_variants.to_csv("acmg_hc_lof_variants.csv", index=False)
    print(f"HC LoF filtering complete. Variants retained: {acmg_lof_variants.shape[0]}")
else:
    print("HC LoF filtering skipped: 'hc_lof' DataFrame was not found.")