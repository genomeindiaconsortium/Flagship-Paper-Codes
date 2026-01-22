import pandas as pd
import glob

# Annotating variants with PharmGKB
merged_df = pd.DataFrame()

snp_info_df = pd.read_excel('pharmgkb_rsids.xlsx')
snp_info_df.drop(columns=['ref', 'alt'], inplace=True)
file_list = glob.glob('pop_freq_alt/*.frq')

for file in file_list:
    df = pd.read_csv(file, sep=r'\s+')
    df = df[['SNP', 'A1', 'A2', 'MAF']]
    df["SNP"] = df["SNP"].str.split(";")
    df = df.explode("SNP")
    maf_col_name = file.replace('\\', '.').split('.')[1]
    df.rename(columns={'MAF': maf_col_name}, inplace=True)

    if merged_df.empty:
        merged_df=df
    else:
        df = df[['SNP', maf_col_name]]
        merged_df=pd.merge(merged_df, df, on='SNP', how='outer')

non_numeric_cols = merged_df.iloc[:, :3]
numeric_cols = merged_df.iloc[:, 3:]
numeric_cols = numeric_cols.reindex(sorted(numeric_cols.columns, key=lambda x:int(x)), axis = 1)
merged_df = pd.concat([non_numeric_cols, numeric_cols], axis=1)

final_df = pd.merge(merged_df,
                    snp_info_df,
                    left_on=merged_df['SNP'],
                    right_on="variant",
                    how="inner")

final_df = final_df[~final_df["chr"].astype(str).str.startswith("H", na=False)]
final_df = final_df.drop_duplicates()
final_df.drop(columns=['variant'], inplace=True)
final_df.rename(columns={'SNP':'variant'}, inplace=True)

cols = ['variant','chr','start','end','A1','A2','gene','type','level of evidence','chemicals','phenotypes'] + [col for col in merged_df.columns if col not in ['SNP', 'A1', 'A2']]
final_df = final_df[cols]

final_df['chr'] = pd.to_numeric(final_df['chr'])
final_df = final_df.sort_values('chr')
final_df.to_excel('variant_allele_freq.xlsx', index=False)

# Calculating genotype counts and frequency of each pharmacogenomic variant for each GI populations
merged_df = pd.DataFrame()
results = []

file_list = glob.glob("het_hom\\*.frqx")
for file in file_list:
    df = pd.read_csv(file, sep=r'\s+')
    df = df[['SNP', 'C(HOM_A1)', 'C(HET)', 'C(HOM_A2)']]
    col_name = file.replace('\\', '.').split('.')[7]
    df = pd.melt(df, id_vars=['SNP'], value_vars=['C(HOM_A1)', 'C(HET)', 'C(HOM_A2)'], var_name="Genotype", value_name=col_name)

    if merged_df.empty:
        merged_df=df
    else:
        merged_df=pd.merge(merged_df, df, on=['SNP', 'Genotype'], how='left')

merged_df.columns = merged_df.columns.astype(str)
melted_df = pd.melt(merged_df, 
                    id_vars=['SNP', 'Genotype'], 
                    value_vars=[str(i) for i in range(1, 83)], 
                    var_name='Populations', 
                    value_name='Sample_Count')

melted_df.to_excel("het_hom_count.xlsx", index=False)

for snp in melted_df['SNP'].unique():
    snp_df = melted_df[melted_df['SNP'] == snp]

    for pop in snp_df['Populations'].unique():
        sub_df = snp_df[snp_df['Populations'] == pop]

        hom_ref = sub_df.loc[sub_df['Genotype'] == 'C(HOM_A2)', 'Sample_Count'].sum()
        het = sub_df.loc[sub_df['Genotype'] == 'C(HET)', 'Sample_Count'].sum()
        hom_alt = sub_df.loc[sub_df['Genotype'] == 'C(HOM_A1)', 'Sample_Count'].sum()

        total_individuals = hom_ref + het + hom_alt

        ref_genotype_freq = hom_ref / total_individuals if total_individuals else 0
        alt_genotype_freq = hom_alt / total_individuals if total_individuals else 0
        het_genotype_freq = het / total_individuals if total_individuals else 0

        results.append({
            'SNP': snp,
            'Population': pop,
            'ref_genotype_freq': ref_genotype_freq,
            'alt_genotype_freq': alt_genotype_freq,
            'het_genotype_freq': het_genotype_freq
        })

freq_df = pd.DataFrame(results)
freq_df.rename(columns={'SNP': 'variant', 'Population':'populations'}, inplace=True)
freq_df.to_excel("het_hom_freq.xlsx", index=False)

# Average contribution of risk-associated genotypes per variant for a given drug within a population
merged = pd.DataFrame()
pivoted = pd.DataFrame()

imp_var = pd.read_excel("imp_var_detailed.xlsx", sheet_name="Imp") # Prioritized list of 49 pharmacogenomic SNPs and indels important for GI
imp_var["populations"] = imp_var["populations"].astype(str)

genotype_count = pd.read_excel("het_hom_count.xlsx")
genotype_count = genotype_count.rename(columns={"SNP": "variant", "Populations": "populations"})
genotype_count["populations"] = genotype_count["populations"].astype(str)
merged = pd.merge(imp_var, genotype_count, on=["variant", "populations"], how="left")
merged.to_excel("imp_var_genotype_count.xlsx", index=False)

pivoted = merged.pivot_table(index=["variant", "imp_allele", "allele", "type", "chemicals", "populations", "Linguistic_Group",
                                    "Admixture_Cluster", "Testing"],
                                    columns="Genotype", values="Sample_Count", fill_value=0).reset_index()
pivoted.columns.name = None
pivoted = pivoted.rename(columns={
    "C(HOM_A1)": "HOM_A1",
    "C(HOM_A2)": "HOM_A2",
    "C(HET)": "HET"
})

def calculate_effect_freq(row):
    total = row['HOM_A1'] + row['HET'] + row['HOM_A2']
    if total == 0:
        return None
    if row['allele'] == 'Ref':
        return (row['HOM_A2'] + row['HET']) / total
    elif row['allele'] == 'Alt':
        return (row['HOM_A1'] + row['HET']) / total
    else:
        return None

pivoted['Effect_Allele_Freq'] = pivoted.apply(calculate_effect_freq, axis=1)
pivoted['Risk_Genotype_Freq'] = pivoted['Effect_Allele_Freq']

chem_pop_freq = pivoted.groupby(['chemicals', 'populations', 'type'])['Effect_Allele_Freq'].sum().reset_index()
variant_counts = pivoted.groupby(['chemicals', 'populations', 'type'])['variant'].nunique().reset_index()
variant_counts = variant_counts.rename(columns={'variant': 'variant_count'})
chem_pop_freq = pd.merge(chem_pop_freq, variant_counts, on=['chemicals', 'populations', 'type'], how='left')
chem_pop_freq['Mean_Risk_Genotype_Freq'] = chem_pop_freq['Effect_Allele_Freq'] / chem_pop_freq['variant_count']
pivoted = pivoted.drop(columns='Effect_Allele_Freq')
chem_pop_freq = chem_pop_freq.drop(columns='Effect_Allele_Freq')
chem_pop_freq = pd.merge(pivoted, chem_pop_freq, on=['chemicals', 'populations', 'type'], how='right')

chem_pop_freq = chem_pop_freq.drop_duplicates()
chem_pop_freq.to_excel("drug_genotype_freq.xlsx", index=False)