import pandas as pd
import numpy as np
from cyvcf2 import VCF
import csv

# Actionable SNPs
risk_map = {}
with open ("11_snp_actionable.txt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        rsid = row["rsid"]
        ref_is_risk = row["ref_is_risk"].lower() == "yes"
        risk_genotype = row["risk_genotype"].lower() == "hom"
        risk_map[rsid] = {"ref_is_risk":ref_is_risk, "risk_genotype":risk_genotype}

vcf = VCF("11_snp_gi.vcf.gz")
samples = vcf.samples
counts = {s:0 for s in samples}

for var in vcf:
    rsid = var.ID
    ref_is_risk = risk_map[rsid]["ref_is_risk"]
    risk_genotype = risk_map[rsid]["risk_genotype"]
    
    for i, gt in enumerate(var.genotypes):
        g1, g2 = gt[0], gt[1]

        if g1 == -1 or g2 == -1:
            continue
        
        if risk_genotype:
            if ref_is_risk: 
                if g1 == 0 and g2 == 0:
                    counts[samples[i]] += 1
            else:
                if g1 == 1 and g2 == 1:
                    counts[samples[i]] += 1
        else:
            if ref_is_risk:
                if g1 == 0 and g2 == 0 or ((g1 == 0 and g2 == 1) or (g1 == 1 and g2 == 0)):
                    counts[samples[i]] += 1
            else:
                if g1 == 1 and g2 == 1 or ((g1 == 0 and g2 == 1) or (g1 == 1 and g2 == 0)):
                    counts[samples[i]] += 1

with open("11_snp_actionable_counts.tsv", "w") as out:
    out.write("Sample\tRisk_SNP_Count\n")
    for sample in samples:
        out.write(f"{sample}\t{counts[sample]}\n")

# VKORC1 genotype calls
merge_df = pd.DataFrame()
final_counts = pd.DataFrame()

risk_snp_counts = pd.read_csv("11_snp_actionable_counts.tsv", sep="\t")
vkorc1_genotype_calls = pd.read_csv("vkorc1_9923231_sample_calls.txt", sep="\t", header=None, names=["Sample", "vkorc1_genotype"])
merge_df = pd.merge(risk_snp_counts, vkorc1_genotype_calls, on="Sample", how="left")

merge_df["Risk_SNP_Count_with_vkorc1"] = np.where(
    (merge_df["vkorc1_genotype"] == "0/1") | (merge_df["vkorc1_genotype"] == "1/1"),
    merge_df["Risk_SNP_Count"] + 1,
    merge_df["Risk_SNP_Count"]
)
merge_df.to_csv("49_snp_actionable_with_vkorc1.tsv", sep="\t", index=False)

# Actionable diplotypes for actionable star alleles
target_alleles = {
    "cyp2b6": ["*6"],
    "cyp2c19": ["*2"],
    "cyp2c9": ["*3", "*11"],
    "tpmt": ["*3A"],
    "ugt1a1" : ["*6"]
}

non_normal_phenotypes = {
    'intermediate_metabolizer',
    'poor_metabolizer',
    'slow_metabolizer',
    'rapid_metabolizer',
    'ultrarapid_metabolizer',
    'increased_function',
    'decreased_function',
    'poor_function'
}

for gene, alleles in target_alleles.items():
    file_path = f"stargazer//output-{gene}.sg-genotype.txt"
    df = pd.read_csv(file_path, sep = "\t", usecols=["name", "hap1_main", "hap2_main", "phenotype"])

    df["has_target"] = df[["hap1_main", "hap2_main"]].apply(lambda row: any(allele in alleles for allele in row), axis=1)
    df["is_non_normal"] = df["phenotype"].isin(non_normal_phenotypes)

    df[f"{gene}_risk"] = df["has_target"] & df["is_non_normal"]

    subset = df[["name", f"{gene}_risk"]].rename(columns={"name": "FID"})
    if final_counts.empty:
        final_counts = subset
    else:
        final_counts = pd.merge(final_counts, subset, on="FID", how="outer")

cyp2d6_phenotypes = pd.read_excel("cyp2d6_analysis\\merged_genotypes.xlsx",
                                  usecols=["Sample", "Genotype", "Metabolizer_Status"])
target_alleles = ["*4", "*2x2"]

non_normal_phenotypes = {
    "Intermediate Metabolizer",
    "Poor Metabolizer",
    "Ultrarapid Metabolizer"
}

cyp2d6_phenotypes[["hap1", "hap2"]] = cyp2d6_phenotypes["Genotype"].str.split("/", expand=True)
cyp2d6_phenotypes["has_target"] = cyp2d6_phenotypes[["hap1", "hap2"]].apply(lambda row: any(allele in target_alleles for allele in row),
                                                                            axis=1)
cyp2d6_phenotypes["is_non_normal"] = cyp2d6_phenotypes["Metabolizer_Status"].isin(non_normal_phenotypes)
cyp2d6_phenotypes["cyp2d6_risk"] = cyp2d6_phenotypes["has_target"] & cyp2d6_phenotypes["is_non_normal"]

final_counts_with_cyp2d6 = pd.DataFrame()

final_counts_with_cyp2d6 = cyp2d6_phenotypes[["Sample", "cyp2d6_risk"]].rename(columns={"Sample":"FID"})
final_counts_with_cyp2d6 = pd.merge(final_counts, final_counts_with_cyp2d6, on="FID", how="outer")

risk_cols = [col for col in final_counts_with_cyp2d6.columns if col.endswith("_risk")]
final_counts_with_cyp2d6[risk_cols] = final_counts_with_cyp2d6[risk_cols].fillna(False)
final_counts_with_cyp2d6["non_normal_phenotype_count"] = final_counts_with_cyp2d6[risk_cols].sum(axis=1)
final_counts_with_cyp2d6.to_csv("non_normal_phenotype_counts_32_star_alleles.tsv", sep="\t", index=False)