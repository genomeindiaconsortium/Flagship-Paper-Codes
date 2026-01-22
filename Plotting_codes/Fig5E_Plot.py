import matplotlib.pyplot as plt
import os
import numpy as np
import time
import pandas as pd
import plotly.express as px
import datatable as dt
import glob


df1 = dt.fread("Results/Metrics_per_variant_SNPS_INDELS.txt")
df = df1.to_pandas()[['ID','Sensitivity','Specificity','MAF','R2','Precision']]
print(df.head())


fig = plt.figure(figsize=(7, 6))
gs  = fig.add_gridspec(nrows=1, ncols=2) 
#gs=fig.add_gridspec(1, 2, width_ratios=[1,1])
axes_main = [fig.add_subplot(gs[0, i]) for i in [0,1]] 
#axes_top = [fig.add_subplot(gs[0, i], sharex=axes_main[i]) for i in [0,1]]

def clean_label(val):
    if val < 0:
        return "0.00"
    else:
        return f"{val:.2f}"

hb=None
#metrics = ['Sensitivity', 'Specificity', 'Precision']
metrics = ['Sensitivity', 'Specificity']
for i, metric in enumerate(metrics):
	#x, y, sizes, colors = generate_plot(metric)
	x=df['R2']
	y=df[metric]
	z=df['MAF']
     
    #Hexbin plot of each metric vs INFO
	hb=axes_main[i].hexbin(x=x, y=y, C=z,reduce_C_function=np.mean, gridsize=20, cmap='cividis', mincnt=1)
	#axes_main[i].set_title(f'{metric} vs INFO')
   
	axes_main[i].set_ylim(0, 1)

	x_bins, x_bin_edges = pd.cut(x, bins=20, retbins=True)
	y_bins, y_bin_edges = pd.cut(y, bins=20, retbins=True)
	
	axes_main[i].tick_params(axis='x')
	axes_main[i].tick_params(left=True, bottom=True, labelleft=True, labelbottom=True,axis='both',labelsize=14)
	axes_main[i].set_xlabel('INFO', fontdict={"fontsize": 16})
	axes_main[i].set_ylabel(f'{metric}',fontdict={"fontsize": 16})
	axes_main[i].set_xlim(0, 1)
	axes_main[i].set_ylim(0, 1)
	# Set ticks at 0.0, 0.1, 0.2, … 1.0
	ticks = np.arange(0, 1.001, 0.2)
	axes_main[i].set_xticks(ticks)
	axes_main[i].set_yticks(ticks)

#Add MAF Colorbar at the bottom
cbar = fig.colorbar(
    hb,
    ax = axes_main,            # share colorbar for both subplots
    orientation = 'horizontal',
    fraction    = 0.05,        # adjusts thickness of the bar
    pad         = 3,         # spacing between plots and colorbar
    label       = 'Mean MAF'
)
cbar.set_label('Mean MAF', fontsize=16)
fig.subplots_adjust(wspace=0.5, hspace=0.5)
plt.tight_layout()
plt.subplots_adjust(right=0.85,bottom =0.3) 
#plt.tight_layout()
# Save the plot to a file
plt.savefig("Results/UKB_ALL_INFO_MAF_plot_MAC3_v1_SNPS_INDELS_28112025.png", bbox_inches='tight',dpi=600)