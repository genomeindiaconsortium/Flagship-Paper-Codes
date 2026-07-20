## README

## Prerequisites & Dependencies
To run the visualization scripts, you will need R and the following packages:
* ggplot2
* ggpubr
* ggrepel
* dplyr

Note: Details regarding the upstream bioinformatics tools used to generate the raw PCA data can be found in Supplementary File S4 of the paper.

## Data Availability
The input data files (pca_2b_data.txt) should be placed in their respective directories (pca_2b/) relative to the scripts. 

## Repository Structure & Code Descriptions

## PCA Plot for Genome India and Neighboring Populations (Figure 2b)
* Script Name: plot_pca_fig2b.R
* Description: Generates Figure 2b, a PCA plot comparing the Genome India population with neighboring populations from the GenomeAsia and HGDP datasets. The visualization uses color to represent linguistic groups or geographic regions, and shape to denote consortium name or tribal status.
* Input Data: pca_2b_data.txt (Expected in a pca_2b/ directory)

