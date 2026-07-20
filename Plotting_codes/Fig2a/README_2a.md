# README

## Prerequisites & Dependencies
To run the visualization scripts, you will need R and the following packages:
* ggplot2
* ggpubr

Note: Details regarding the upstream bioinformatics tools used to generate the raw PCA data can be found in Supplementary File S4 of the paper.

## Data Availability
The input data files (e.g., pca_2a_data.csv) should be placed in their respective directories (e.g., pca_2a/) relative to the scripts. 

## Repository Structure & Code Descriptions
## PCA Plot for Population Structure (Figure 2a)
* Script Name: plot_pca_fig2a.R
* Description: Generates Figure 2a, a Principal Component Analysis (PCA) plot visualizing the population structure of all samples. The plot maps individuals based on their linguistic affiliation (represented by color combinations) and tribal status (represented by shape).
* Input Data: pca_2a_data.csv (Expected in a pca_2a/ directory)

