import pandas as pd

# Annotating locations to variants of pharmGKB
pharmgkb = pd.read_excel("clinicalVariants.xlsx")
pharmgkb_rsids = pharmgkb.loc[pharmgkb['variant'].str.startswith('rs')]
pharmgkb_loc = pd.read_csv("vep_result.txt", sep="\t", usecols=['#Uploaded_variation',
                                                                  'Location',
                                                                  'Allele',
                                                                  'REF_ALLELE'])

df = []

df = pd.merge(pharmgkb_rsids,
              pharmgkb_loc,
              left_on='variant',
              right_on='#Uploaded_variation',
              how="left")

df = df.rename(columns={'Allele':'alt', 'REF_ALLELE':'ref'})
df[['chr','start','end']] = df['Location'].str.split('[:-]', expand=True)
df[['chr','start','end']] = df[['chr','start','end']].astype(str)
cols = ['chr', 'start', 'end', 'ref', 'alt', 'variant', 'gene', 'type', 'level of evidence', 
        'chemicals', 'phenotypes']
df = df[cols]
df = df.drop_duplicates()
df.to_excel('pharmgkb_rsids.xlsx', index=False)

# Fetch star alleles

pharmgkb_star_alleles = pharmgkb.loc[~pharmgkb['variant'].str.startswith('rs')]
df = pd.DataFrame(pharmgkb_star_alleles)
df['variant'] = df['variant'].str.split(',')
df = df.explode('variant')
df.to_excel("pharmgkb_haplotypes.xlsx", index=False)

# Fetching testing levels and their prescribing information
drug_annotations = pd.read_csv(
    "drugLabels.tsv", 
    usecols=[
    "Testing Level", 
    "Has Prescribing Info", 
    "Chemicals", 
    "Genes", 
    "Variants/Haplotypes"],
sep=r'\t').drop_duplicates()

drug_annotations["Variants/Haplotypes"] = drug_annotations["Variants/Haplotypes"].str.split("; ")
drug_annotations = drug_annotations.explode("Variants/Haplotypes").reset_index(drop=True)
drug_annotations = drug_annotations[
    drug_annotations["Variants/Haplotypes"].notna() &
    drug_annotations["Variants/Haplotypes"].str.strip().ne("") &
    drug_annotations["Variants/Haplotypes"].str.strip().ne("None")
]
drug_annotations = drug_annotations[drug_annotations["Testing Level"].isin([
    "Testing Required", "Testing Recommended", "Actionable PGx"])]
drug_annotations = drug_annotations.sort_values("Testing Level")
drug_annotations.to_excel("imp_testing_levels.xlsx", index=False)