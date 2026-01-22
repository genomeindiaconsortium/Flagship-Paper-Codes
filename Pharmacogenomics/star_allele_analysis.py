import pandas as pd
import numpy as np
import glob
import os

populations = pd.read_excel("pop_info.xlsx", usecols=["FID", "Population", "Population_Order"])

# Stargazer analysis
# Output from Stargazer
genotype_files = glob.glob("stargazer\\output-*.sg-genotype.txt")

# Star-allele frequency of GI and their comparison to 1KGP3 populations
results = []
def calculate_allele_frequency(df, hap1_col, hap2_col):
    hap_concat=pd.concat([df[hap1_col],df[hap2_col]],axis=0)
    allele_counts = hap_concat.value_counts()
    total_alleles = allele_counts.sum()
    allele_freq = allele_counts / total_alleles
    return allele_freq

thousand_genomes_df = pd.read_csv("1000_Genome_Pharmacogenomics.csv")

def explode_haplotypes(df, hap1_col, hap2_col):
    df[hap1_col] = df[hap1_col].str.rstrip(";").str.split(";")
    df[hap2_col] = df[hap2_col].str.rstrip(";").str.split(";")
    df = df.explode(hap1_col).explode(hap2_col)
    return df

thousand_genomes_df = explode_haplotypes(thousand_genomes_df, 'Haplotype1', 'Haplotype2')

for file in genotype_files:
    gene_name = file.split("stargazer\\output-")[1].split(".sg-genotype.txt")[0]
    df = pd.read_csv(file, sep="\t")
    df = df.merge(populations, left_on='name', right_on='FID', how='left')
    df = df[~(df["Population"] == "Siddi")]

    allele_freqs = calculate_allele_frequency(df, 'hap1_main', 'hap2_main')
    indian_allele_freq = allele_freqs.to_dict()

    results.append({"Gene": gene_name, "Genome India": indian_allele_freq})

    gene_data = thousand_genomes_df[thousand_genomes_df['Gene'] == gene_name.upper()]
    pop_filtered_allele_freqs = {"ALL": {}}

    genome_allele_freqs = calculate_allele_frequency(gene_data, 'Haplotype1', 'Haplotype2')
    filtered_genome_allele_freqs = {allele: genome_allele_freqs[allele] for allele in indian_allele_freq.keys() if allele in genome_allele_freqs}
    pop_filtered_allele_freqs["ALL"] = filtered_genome_allele_freqs

    superpopulations = ["AFR", "AMR", "EUR", "EAS", "SAS"]
    for pop in superpopulations:
        pop_data = gene_data[gene_data['Superpopulation code'] == pop]
        pop_allele_freqs = calculate_allele_frequency(pop_data, 'Haplotype1', 'Haplotype2')
        filtered_pop_allele_freqs = {allele: pop_allele_freqs[allele] for allele in indian_allele_freq.keys() if allele in pop_allele_freqs}
        pop_filtered_allele_freqs[pop] = filtered_pop_allele_freqs
    
    results.append({"Gene": gene_name, "Genome India": indian_allele_freq, "1000G Allele Freq": pop_filtered_allele_freqs})

table_results = []
for entry in results:
    gene_name = entry["Gene"]
    for allele in entry["Genome India"].keys():
        revised_gene_name = f"{gene_name}{allele}"
        row = {"Gene": revised_gene_name, "Genome India": entry["Genome India"].get(allele, 0)}
        
        if "1000G Allele Freq" in entry:
            for pop in ["ALL", "AFR", "AMR", "EUR", "EAS", "SAS"]:
                row[f"1KGP3 {pop}"] = entry["1000G Allele Freq"].get(pop, {}).get(allele, 0)
        else:
            for pop in ["ALL", "AFR", "AMR", "EUR", "EAS", "SAS"]:
                row[f"1KGP3 {pop}"] = 0
        
        if row["Genome India"] > 0 and row["1KGP3 ALL"] > 0:
            table_results.append(row)

results_df = pd.DataFrame(table_results)
results_df.to_excel("allele_frequencies.xlsx", index=False)

# Population-specific star-allele frequency of GI
merged_data = populations.copy()
gene_names = []
final_data = pd.DataFrame()

for file in genotype_files:
    gene_name = file.split("stargazer\\output-")[1].split(".sg-genotype.txt")[0]
    gene_names.append(gene_name)

    df = pd.read_csv(file, sep="\t", usecols=["name", "hap1_main", "hap2_main"])
    df.rename(columns={"hap1_main": f"{gene_name}_hap1", "hap2_main": f"{gene_name}_hap2", "name" : "FID"}, inplace=True)

    merged_data = pd.merge(merged_data,
                           df,
                           on="FID",
                           how="left")

merged_data = merged_data[~(merged_data["Population"] == "Siddi")]
merged_data = merged_data.drop_duplicates()
pop_id = merged_data["Population_Order"].unique()
merged_data.to_excel("merged_genotypes.xlsx", index=False)

for pop in pop_id:
    filtered_pop_df = merged_data.loc[merged_data["Population_Order"] == pop]
    for gene in gene_names:
        hap1 = f"{gene}_hap1"
        hap2 = f"{gene}_hap2"

        if hap1 not in filtered_pop_df.columns or hap2 not in filtered_pop_df.columns:
            continue

        allele_freqs = calculate_allele_frequency(filtered_pop_df, hap1, hap2)

        for allele, frequency in allele_freqs.items():
            gene_allele = f"{gene}{allele}"
            final_data.loc[gene_allele, pop] = frequency

final_data.reset_index(inplace=True)
final_data.rename(columns={"index": "Star_Allele"}, inplace=True)
final_data = final_data.drop_duplicates()
final_data = final_data.fillna(0)
final_data.to_excel("pop_specific_allele_frequencies.xlsx", index=False)

# Prioritizing important star-alleles
indian_global = pd.DataFrame()
final_merged = pd.DataFrame()

global_freq = pd.read_excel("allele_frequencies.xlsx", usecols=["Gene", "1KGP3 SAS"]).rename(columns={"1KGP3 SAS":"SAS", 
                                                                                                      "Gene":"Star_Allele"})
global_freq["Star_Allele"] = global_freq["Star_Allele"].str.upper()

pop_spec_freq = pd.read_excel("pop_specific_allele_frequencies.xlsx").melt(id_vars="Star_Allele",
                                                                           value_vars=[i for i in range(1,83)],
                                                                           var_name="Populations",
                                                                           value_name="Frequency")
pop_spec_freq["Star_Allele"] = pop_spec_freq["Star_Allele"].str.upper()

haplotype_info = pd.read_excel("pharmgkb_haplotypes.xlsx")[lambda haplotype_info:haplotype_info["level of evidence"] != 4].drop(
    ["phenotypes", "gene", "level of evidence"], axis=1
    ).rename(columns = {"variant":"Star_Allele"})
haplotype_info["Star_Allele"] = haplotype_info["Star_Allele"].str.strip().str.upper()
indian_global = pd.merge(global_freq, pop_spec_freq, on="Star_Allele", how="inner")
final_merged = pd.merge(indian_global, haplotype_info, on="Star_Allele", how="left")
final_merged["compare"] = final_merged["Frequency"] - final_merged["SAS"]
final_merged.to_excel("imp_star_alleles.xlsx", index=False)

# Cyrius analysis for CYP2D6 star-alleles
# Assembling output of Cyrius from four nodal centres
df = pd.DataFrame()
cyrius_result = pd.DataFrame()

input_directory = "igib" # cbr, ccmb, igib, nibmg
output_file = "cyp2d6_igib.xlsx"

# CBR
cyrius_result = pd.read_csv(input_directory, sep="\t", keep_default_na=False)[lambda cyrius_result:cyrius_result["Filter"] == "PASS"]

# CCMB, IGIB
for files in os.listdir(input_directory):
    if files.endswith("tsv"):
        result_path = os.path.join(input_directory, files)

        df = pd.read_csv(result_path, sep="\t", keep_default_na=False)[lambda cyrius_result:cyrius_result["Filter"] == "PASS"]
        if cyrius_result.empty:
            cyrius_result = df
        else:
            cyrius_result = pd.concat([cyrius_result, df], ignore_index=True)

cyrius_result.drop_duplicates()
cyrius_result.to_excel(output_file, index=False)

# CYP2D6 frequency calculation
cbr = pd.read_excel("cyp2d6_cbr.xlsx", usecols=["Sample", "Genotype"])
ccmb = pd.read_excel("cyp2d6_ccmb.xlsx", usecols=["Sample", "Genotype"])
igib = pd.read_excel("cyp2d6_igib.xlsx", usecols=["Sample", "Genotype"])
nibmg = pd.read_excel("cyp2d6_nibmg.xlsx", usecols=["Sample", "Genotype"])
gi = pd.read_excel("pop_info.xlsx", dtype="str")[lambda gi:gi["Population"] != "Siddi"]

cyp2d6_gi = pd.concat([cbr, ccmb, igib, nibmg])
cyp2d6_gi = cyp2d6_gi.drop_duplicates()

cyp2d6_phenotype = pd.read_excel("cyp2d6_phenotype/CYP2D6_Diplotype_Phenotype_Table.xlsx",
                                 usecols=["CYP2D6 Diplotype", "Coded Diplotype/Phenotype Summary"])
cyp2d6_phenotype["Coded Diplotype/Phenotype Summary"] = cyp2d6_phenotype["Coded Diplotype/Phenotype Summary"].str.replace("CYP2D6", "", regex=False).str.strip()

cyp2d6_phenotype["reversed_diplotype"] = cyp2d6_phenotype["CYP2D6 Diplotype"].str.split("/").apply(lambda x:f"{x[1]}/{x[0]}")
cyp2d6_phenotype_final = cyp2d6_phenotype.melt(id_vars="Coded Diplotype/Phenotype Summary",
                                               value_vars=["CYP2D6 Diplotype", "reversed_diplotype"],
                                               value_name="Diplotype").drop("variable", axis=1)

merged = pd.merge(gi, cyp2d6_gi, left_on="FID", right_on="Sample", how="inner")
merged = merged.drop_duplicates()
merged_phenotypes = pd.merge(merged, cyp2d6_phenotype_final, left_on="Genotype", right_on="Diplotype", how="left")
merged_phenotypes = merged_phenotypes.drop_duplicates()
merged_phenotypes["Coded Diplotype/Phenotype Summary"] = merged_phenotypes["Coded Diplotype/Phenotype Summary"].replace(np.NaN, "Unknown Metabolizer")
merged_phenotypes = merged_phenotypes.rename(columns={"Coded Diplotype/Phenotype Summary":"Metabolizer_Status"})
merged_phenotypes = merged_phenotypes.drop("Diplotype", axis=1)
merged_phenotypes.to_excel("merged_genotypes.xlsx", index=False)
merged_phenotypes[["hap1", "hap2"]] = merged_phenotypes["Genotype"].str.split("/", expand=True)

# Overall CYP2D6 star-allele frequency in GI
df = pd.read_excel("merged_genotypes.xlsx", usecols=["Sample", "Genotype"])
df[["hap1", "hap2"]] = df["Genotype"].str.split("/", expand=True)
overall_allele_freq = calculate_allele_frequency(df, "hap1", "hap2")
rows = []

for allele, freq in overall_allele_freq.items():
    rows.append({"cyp2d6_allele":allele, "cyp2d6_freq":freq})

overall_freq_df = pd.DataFrame(rows)
overall_freq_df.to_excel("overall_cyp2d6_freq.xlsx", index=False)

# Population-specific CYP2D6 star-allele frequency in GI
rows = []

for pop in merged_phenotypes["Population_Order"].unique():
    filtered_df = merged_phenotypes[merged_phenotypes["Population_Order"] == pop]
    allele_freq = calculate_allele_frequency(filtered_df, "hap1", "hap2")

    for allele, freq in allele_freq.items():
        rows.append({"Population_Order":pop, "cyp2d6_allele":allele, "cyp2d6_freq":freq})

cyp2d6_freq = pd.DataFrame(rows)
final_df = pd.merge(gi, cyp2d6_freq, on="Population_Order", how="inner")
final_df = final_df.drop(columns=["FID", "IID"])
final_df["Population_Order"] = final_df["Population_Order"].astype(int)
final_df = final_df.drop_duplicates()

pivot = final_df.pivot_table(index='cyp2d6_allele', 
                            columns='Population_Order',
                            values='cyp2d6_freq',
                            aggfunc='sum',
                            fill_value=0)

pivot = pivot.reindex(columns=range(1, 83), fill_value=0)
pivot = pivot.reset_index()
pivot.to_excel('pop_spec_cyp2d6_freq.xlsx', index=False)