import os
import numpy as np
import time
import pandas as pd
import plotly.express as px
import datatable as dt
import glob


conf_paths=glob.glob("Results/Sensitivity_7628_9768_ConfusionMetrics_*.csv_reheader.csv")


conf_list = []
for file_path in conf_paths:
    df = dt.fread(file_path)
    df1 = df.to_pandas()
    conf_list.append(df1)

#print(df_list)
data = pd.concat(conf_list, ignore_index=True)

print("Read and Concatenated!")
data.to_csv("Results/Sensitivity_7628_9768_ConfusionMetrics_1-22.csv",index=False)


info_paths=glob.glob("../OUTPUT/*/UKB_IND_PAK_BAN-*-GI-unfiltered-MAF-MAC-Metrics.txt.gz")

df_list = []
for file_path in info_paths:
    df = dt.fread(file_path,header=None)
    df1 = df.to_pandas()
    df_list.append(df1)

#print(df_list)
INFO = pd.concat(df_list, ignore_index=True)
INFO.columns=["ID","MAF","MAC","R2"]


#INFO=pd.read_csv("22/SANSCOG-2_108-22-GI-unfiltered-MAF-MAC-Metrics.txt",delim_whitespace=True)
#print(INFO.head())
Conf=data.merge(INFO, on="ID", how="left")
print(Conf.head())

#Conf=data
print("Merged INFO.")
print("Calculating metrics row-wise....")
Conf['Sensitivity_HomR']=Conf['TP_HomR']/(Conf['TP_HomR']+Conf['FN_HomR'])
Conf['Sensitivity_Het']=Conf['TP_Het']/(Conf['TP_Het']+Conf['FN_Het'])
Conf['Sensitivity_HomA']=Conf['TP_HomA']/(Conf['TP_HomA']+Conf['FN_HomA'])
 
Conf['Specificity_HomR']=Conf['TN_HomR']/(Conf['TN_HomR']+Conf['FP_HomR'])
Conf['Specificity_Het']=Conf['TN_Het']/(Conf['TN_Het']+Conf['FP_Het'])
Conf['Specificity_HomA']=Conf['TN_HomA']/(Conf['TN_HomA']+Conf['FP_HomA'])

Conf['Precision_HomR']=Conf['TP_HomR']/(Conf['TP_HomR']+Conf['FP_HomR'])
Conf['Precision_Het']=Conf['TP_Het']/(Conf['TP_Het']+Conf['FP_Het'])
Conf['Precision_HomA']=Conf['TP_HomA']/(Conf['TP_HomA']+Conf['FP_HomA'])


#MicroAveraging
Conf['Sensitivity']=(Conf['TP_HomR']+Conf['TP_Het']+Conf['TP_HomA'])/((Conf['TP_HomR']+Conf['TP_Het']+Conf['TP_HomA'])+(Conf['FN_HomR']+Conf['FN_Het']+Conf['FN_HomA']))
Conf['Specificity']=(Conf['TN_HomR']+Conf['TN_Het']+Conf['TN_HomA'])/((Conf['TN_HomR']+Conf['TN_Het']+Conf['TN_HomA'])+(Conf['FP_HomR']+Conf['FP_Het']+Conf['FP_HomA']))
Conf['Precision']=(Conf['TP_HomR']+Conf['TP_Het']+Conf['TP_HomA'])/((Conf['TP_HomR']+Conf['TP_Het']+Conf['TP_HomA'])+(Conf['FP_HomR']+Conf['FP_Het']+Conf['FP_HomA']))
print("Done.")


Conf.to_csv("Results/Metrics_per_variant_SNPS_INDELS.txt",index=False,sep="\t")
'''
prop=len(Conf[(Conf['Sensitivity']>0.7) & (Conf['Specificity']>0.7) & (Conf['R2']>0.7)])/len(Conf)
got=len(Conf[(Conf['Sensitivity']>0.7) & (Conf['Specificity']>0.7) & (Conf['R2']>0.7) & (Conf['MAF']<0.01)])/len(Conf[(Conf['Sensitivity']>0.7) & (Conf['Specificity']>0.7) & (Conf['R2']>0.7)])
got2=len(Conf[(Conf['Sensitivity']>0.7) & (Conf['Specificity']>0.7) & (Conf['R2']>0.7) & (Conf['MAF']>0.05)])/len(Conf[(Conf['Sensitivity']>0.7) & (Conf['Specificity']>0.7) & (Conf['R2']>0.7)])
 
'''
print("Calculating overall metrics....")
#Overall Metrics
TP_HomR = Conf['TP_HomR'].sum()
FP_HomR = Conf['FP_HomR'].sum()
TN_HomR = Conf['TN_HomR'].sum()
FN_HomR = Conf['FN_HomR'].sum()
Sensitivity_HomR = TP_HomR/(TP_HomR+FN_HomR)
Specificity_HomR = TN_HomR/(TN_HomR+FP_HomR)
Precision_HomR = TP_HomR/(TP_HomR+FP_HomR)

TP_Het=Conf['TP_Het'].sum()
FP_Het=Conf['FP_Het'].sum()
TN_Het=Conf['TN_Het'].sum()
FN_Het=Conf['FN_Het'].sum()
Sensitivity_Het = TP_Het/(TP_Het+FN_Het)
Specificity_Het = TN_Het/(TN_Het+FP_Het)
Precision_Het = TP_Het/(TP_Het+FP_Het)

TP_HomA=Conf['TP_HomA'].sum()
FP_HomA=Conf['FP_HomA'].sum()
TN_HomA=Conf['TN_HomA'].sum()
FN_HomA=Conf['FN_HomA'].sum()
Sensitivity_HomA = TP_HomA/(TP_HomA+FN_HomA)
Specificity_HomA = TN_HomA/(TN_HomA+FP_HomA)
Precision_HomA = TP_HomA/(TP_HomA+FP_HomA)

#MicroAveraging
Overall_Sensitivity=(TP_HomR+TP_Het+TP_HomA)/(TP_HomR+TP_Het+TP_HomA + FN_HomR+FN_Het+FN_HomA)
Overall_Specificity=(TN_HomR+TN_Het+TN_HomA)/(TN_HomR+TN_Het+TN_HomA + FP_HomR+FP_Het+FP_HomA)
Overall_Precision=(TP_HomR+TP_Het+TP_HomA)/(TP_HomR+TP_Het+TP_HomA + FP_HomR+FP_Het+FP_HomA)

print("Done")

print("Writing overall metrics as output...")

with open('Results/UKB_7628_9768_Overall_Metrics_21_22_SNPS_INDELS.txt', 'w') as f:
	print("\nSensitivity for 0/0:",Sensitivity_HomR, file=f)
	print("\nSpecificity for 0/0:",Specificity_HomR, file=f)
	print("\nPrecision for 0/0:",Precision_HomR, file=f)	

	print("\nSensitivity for 0/1:",Sensitivity_Het, file=f)
	print("\nSpecificity for 0/1:",Specificity_Het, file=f)
	print("\nPrecision for 0/1:",Precision_Het, file=f)  

	print("\nSensitivity for 1/1:",Sensitivity_HomA, file=f)
	print("\nSpecificity for 1/1:",Specificity_HomA, file=f)
	print("\nPrecision for 1/1:",Precision_HomA, file=f)  

	print("\nOverall Sensitivity:",Overall_Sensitivity, file=f)
	print("\nOverall Specificity:",Overall_Specificity, file=f)
	print("\nOverall Precision:",Overall_Precision, file=f)  
	
print("Done")

