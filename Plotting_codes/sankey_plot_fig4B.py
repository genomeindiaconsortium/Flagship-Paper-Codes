import pandas as pd
import numpy as np
import plotly.graph_objects as go

# === Sankey Setup ===
df_sankey= pd.read_table(‘sankey_df.txt’)
“””
###sankey_df.txt looks like this
Consequence	Gene	Type	Phenotype
stop_gained	TTN	LoF	Cardiomyopathy
splice_acceptor	KCNQ1	LoF	Long QT syndrome
stop_gained	MYBPC3	LoF	Cardiomyopathy
stop_gained	PKP2	LoF	Cardiomyopathy
splice_acceptor	ATP7B	LoF	Wilson disease
stop_gained	ATP7B	LoF	Wilson disease
stop_gained	DSC2	LoF	Arrhythmogenic right ventricular dysplasia
stop_gained	LDLR	LoF	Hypercholesterolemia
stop_gained	LDLR	LoF	Hypercholesterolemia
splice_donor	LDLR	LoF	Hypercholesterolemia
splice_donor	RPE65	LoF	Retinal dystrophy

“””


df_sankey_filtered = df_sankey.copy()
if "Variant_Count" not in df_sankey_filtered.columns:
    df_sankey_filtered["Variant_Count"] = 1

sankey_links = df_sankey_filtered.groupby(["Consequence", "Type", "Gene", "Phenotype"]).size().reset_index(name="Count")

conseqs = df_sankey_filtered["Consequence"].unique().tolist()
types = sankey_links["Type"].unique().tolist()
genes = sankey_links["Gene"].unique().tolist()
diseases = sankey_links["Phenotype"].unique().tolist()
dummy_node_label = " "

all_nodes = conseqs + types + genes + diseases + [dummy_node_label]

def italicize(label):
    return f"<i>{label}</i>" if label in genes else label

display_labels = [italicize(label) for label in all_nodes]
label_to_idx = {label: idx for idx, label in enumerate(all_nodes)}

# Category Colors
category_colors = {
    "Eye Disorders": "#9C009C",
    "Cardiac Disorders": "#d62728",
    "Metabolic Disorders": "#2ca02c",
    "Cancer Predisposition Syndromes": "#ffdf00",
    "Genetic Connective Tissue Disorders": "#ff7f0e"
}

gene_category_map = {
    "RPE65": "Eye Disorders",
    "MYBPC3": "Cardiac Disorders", "TTN": "Cardiac Disorders", "PKP2": "Cardiac Disorders",
    "KCNQ1": "Cardiac Disorders", "LDLR": "Cardiac Disorders", "APOB": "Cardiac Disorders",
    "ATP7B": "Metabolic Disorders", "GAA": "Metabolic Disorders", "HFE": "Metabolic Disorders", "BTD": "Metabolic Disorders",
    "BRCA1": "Cancer Predisposition Syndromes", "BRCA2": "Cancer Predisposition Syndromes", "PALB2": "Cancer Predisposition Syndromes",
    "COL3A1": "Genetic Connective Tissue Disorders","DSC2" : "Cardiac Disorders"
}

disease_category_map = {
    "Retinitis pigmentosa": "Eye Disorders",
    "Retinal dystrophy": "Eye Disorders",
    "Cardiomyopathy": "Cardiac Disorders",
    "Long QT syndrome": "Cardiac Disorders",
    "Hypercholesterolemia": "Cardiac Disorders",
    "Wilson disease": "Metabolic Disorders",
    "Glycogen storage disease": "Metabolic Disorders",
    "Hemochromatosis": "Metabolic Disorders",
    "Biotinidase deficiency": "Metabolic Disorders",
    "Breast cancer": "Cancer Predisposition Syndromes",
    "Ehlers-Danlos syndrome": "Genetic Connective Tissue Disorders",
    "Arrhythmogenic right ventricular dysplasia" : "Cardiac Disorders"
}

# Color nodes by category
node_colors = []
for label in all_nodes:
    if label in gene_category_map:
        node_colors.append(category_colors[gene_category_map[label]])
    elif label in disease_category_map:
        node_colors.append(category_colors[disease_category_map[label]])
    elif label == dummy_node_label:
        node_colors.append("rgba(0,0,0,0)")
    else:
        node_colors.append("lightgray")

# Sankey links
sources, targets, values, link_color_array = [], [], [], []

# Consequence → Type
df_counts = (
    df_sankey_filtered.groupby(["Consequence", "Type"], as_index=False)["Variant_Count"].sum()
)

for _, row in df_counts.iterrows():
    cons, typ, vc = row["Consequence"], row["Type"], int(row["Variant_Count"])
    if cons not in label_to_idx or typ not in label_to_idx:
        continue
    sources.append(label_to_idx[cons])
    targets.append(label_to_idx[typ])
    values.append(vc)
    link_color_array.append("rgba(200, 200, 200, 0.6)")

# Type → Gene
for _, row in df_sankey_filtered.groupby(["Type", "Gene"]).size().reset_index(name="Count").iterrows():
    sources.append(label_to_idx[row["Type"]])
    targets.append(label_to_idx[row["Gene"]])
    values.append(row["Count"])
    link_color_array.append('rgba(200, 200, 200, 0.6)')

# Gene → Phenotype
for _, row in df_sankey_filtered.groupby(["Gene", "Phenotype"]).size().reset_index(name="Count").iterrows():
    sources.append(label_to_idx[row["Gene"]])
    targets.append(label_to_idx[row["Phenotype"]])
    values.append(row["Count"])
    link_color_array.append('rgba(200, 200, 200, 0.6)')

# Phenotype → Dummy (to keep nodes aligned)
dummy_idx = label_to_idx[dummy_node_label]
for disease in diseases:
    sources.append(label_to_idx[disease])
    targets.append(dummy_idx)
    values.append(1e-9)
    link_color_array.append('rgba(0,0,0,0)')

# Build Sankey Figure
fig = go.Figure(data=[go.Sankey(
    arrangement="snap",
    node=dict(
        pad=20,
        thickness=20,
        line=dict(color="black", width=0.5),
        label=display_labels,
        color=node_colors,
        hovertemplate="%{label}"
    ),
    link=dict(
        source=sources,
        target=targets,
        value=values,
        color=link_color_array
    )
)])

fig.update_layout(
    height=900,
    width=2100,
    title="<b>Variant Functional Impact and Medical Relevance</b>",
    title_x=0.5,
    font=dict(family="Arial", size=28),
    template="plotly_white",
    margin=dict(t=100, b=70, l=50, r=50),
    showlegend=False
)
fig.write_image("Sankey_final_5.png", scale=4)
fig.show()