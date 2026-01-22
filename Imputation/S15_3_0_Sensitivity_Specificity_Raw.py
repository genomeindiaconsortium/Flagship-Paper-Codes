import os
import pandas as pd
import gzip as gz
import numpy as np
from mpi4py import MPI
import time
import multiprocessing as mp
from tqdm import tqdm
import pysam
import sys

#================================================================
#Function to read vcf file removing lines starting with '##'
#================================================================

def read_vcf(path):
    vcf_in = pysam.VariantFile(path)
    sample_names = list(vcf_in.header.samples)
    n_samp= len(sample_names)
    records = []
    for rec in vcf_in.fetch():
        record = {
            'CHROM': rec.chrom,
            'POS': rec.pos,
            'ID': rec.id,
            'REF': rec.ref,
            'ALT': ','.join(str(a) for a in rec.alts),
            'QUAL': rec.qual,
            'FILTER': ';'.join(rec.filter.keys()) if rec.filter else None,
            'INFO': str(rec.info)
        }

        for sample in sample_names:
            # Preserve the actual GT string including phasing
            gt_str = rec.samples[sample].get('GT')
            phased = rec.samples[sample].phased

            if gt_str is not None:
                sep = '|' if phased else '/'
                # Handle missing alleles represented as None
                allele_strs = [str(a) if a is not None else '.' for a in gt_str]
                record[f'{sample}_GT'] = sep.join(allele_strs)
            else:
                record[f'{sample}_GT'] = './.'

        records.append(record)

    return pd.DataFrame(records), n_samp


#=========================================================================================================================================================================================
#Function to get 3x3 confusion matrix:
#
#   	          		Imputed(a)
#		 	0|0     0|1 or 1|0      1|1
#   --------------------------------------------------
#Actual(b)	0/0     0,0     0,1             0,2     TP=(0,0),       FN=(0,1)+(0,2), FP=(1,0)+(2,0), TN=(1,1)+(1,2)+(2,1)+(2,2),     Sensitivity=TP/TP+FN    Specificity=TN/TN+FP
#      		0/1     1,0     1,1             1,2     TP=(1,1),       FN=(1,0)+(1,2), FP=(0,1)+(2,1), TN=(0,0)+(0,2)+(2,0)+(2,2)      
#   		1/1     2,0     2,1             2,2     TP=(2,2),       FN=(2,0)+(2,1), FP=(0,2)+(1,2), TN=(0,0)+(0,1)+(1,0)+(1,1)      
#

#=======================================================================================================================================================================================

def get_metrics(a,b,n_samp):
	frame= {'Imputed genotypes': a[-n_samp:],'WGS genotypes': b[-n_samp:]}
	df = pd.DataFrame(frame)
	#print(df)

	# Define sets for clarity
	HOM_REF = ['0/0', '0|0']
	HET = ['0/1', '0|1', '1|0']
	HOM_ALT = ['1/1', '1|1']

	
	# Precompute boolean masks
	imp = df['Imputed genotypes']
	wgs = df['WGS genotypes']
	
	# Confusion matrix 3x3:Initialize
	mat=[[0,0,0],[0,0,0],[0,0,0]]

	mat[0][0] = ((imp == '0|0') & wgs.isin(HOM_REF)).sum()
	mat[1][1] = ((imp.isin(['0|1', '1|0'])) & wgs.isin(HET)).sum()
	mat[2][2] = ((imp == '1|1') & wgs.isin(HOM_ALT)).sum()

	mat[0][1] = ((imp.isin(['0|1', '1|0'])) & wgs.isin(HOM_REF)).sum()
	mat[0][2] = ((imp == '1|1') & wgs.isin(HOM_REF)).sum()

	mat[1][0] = ((imp == '0|0') & wgs.isin(HET)).sum()
	mat[1][2] = ((imp == '1|1') & wgs.isin(HET)).sum()

	mat[2][0] = ((imp == '0|0') & wgs.isin(HOM_ALT)).sum()
	mat[2][1] = ((imp.isin(['0|1', '1|0'])) & wgs.isin(HOM_ALT)).sum()
        #print(mat)       
#------------------------Confusion metrics----------------
	TP_HomR=mat[0][0]
	TP_Het=mat[1][1]
	TP_HomA=mat[2][2]
#
	FN_HomR=mat[0][1]+mat[0][2]
	FN_Het=mat[1][0]+mat[1][2]
	FN_HomA=mat[2][0]+mat[2][1]
#
	FP_HomR=mat[1][0]+mat[2][0]
	FP_Het=mat[0][1]+mat[2][1]
	FP_HomA=mat[0][2]+mat[1][2]
#
	TN_HomR=mat[1][1]+mat[1][2]+mat[2][1]+mat[2][2]
	TN_Het=mat[0][0]+mat[0][2]+mat[2][0]+mat[2][2]
	TN_HomA=mat[0][0]+mat[0][1]+mat[1][0]+mat[1][1]
#
	#concat Variants
	new = pd.concat([a[:5], pd.Series([TP_HomR, TP_Het, TP_HomA,FN_HomR, FN_Het, FN_HomA,FP_HomR, FP_Het, FP_HomA,TN_HomR, TN_Het, TN_HomA])],ignore_index=True)
	return new
#End of function



#==========================================================
# Function to chunk sites.txt into lists of sites (chr, pos)
#========================================================
def chunk_sites(input_file, chunk_size=10000):
    sites = []
    with open(input_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Assuming sites.txt format: "chr:pos:ref:alt"
            parts = line.split(':')
            chrom, pos = parts[0], int(parts[1])
            sites.append((chrom, pos))
    for i in range(0, len(sites), chunk_size):
        yield sites[i:i+chunk_size]


def parallel_process_df(imputed_df, wgs_df, n_chunks=None):
	if n_chunks is None:
		n_chunks = cpu_count()

	combined_df = pd.DataFrame({'Imputed genotypes': imputed_df,'WGS genotypes': wgs_df})

	chunk_size = int(np.ceil(len(combined_df) / n_chunks))
	chunks = [combined_df.iloc[i*chunk_size:(i+1)*chunk_size] for i in range(n_chunks)]

	from multiprocessing import Pool
	with Pool(n_chunks) as pool:
		results = pool.map(compute_metrics_chunk, chunks)

	# Sum metrics across chunks
	final_metrics = pd.DataFrame(results).sum()
	return final_metrics

#end of parallel function

#================================================================================
#Main starts from here
#=================================================================================

start=time.time()
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

#chrom=["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"]

# Parse --chrom from command line
if "--chrom" in sys.argv:
    chrom = sys.argv[sys.argv.index("--chrom") + 1]
print(f"Running for chromosome {chrom}")

chunk=range(0,62)

# Get starting offset if running batches
start_rank = int(sys.argv[sys.argv.index("--start-rank") + 1]) if "--start-rank" in sys.argv else 0
global_rank = start_rank + rank

#7628 samples chunked into 62 batches and then each batch was processed separately and then integrated 
#with sum of the metrices per variant.

TOTAL_CHUNKS = 62
if global_rank >= TOTAL_CHUNKS:
    sys.exit(0)

chunk_id = f"{int(chunk[global_rank]):02d}"
#Read datasets imputed and WGS
overlap_imputed,n_samp=read_vcf(f"Chunks/"+chrom+"_IMP_chunk_"+str(chunk_id)+".vcf.gz")
#WGS has genotyped data
overlap_WGS,n_samp=read_vcf(f"Chunks/"+chrom+"_WGS_chunk_"+str(chunk_id)+".vcf.gz")

print(overlap_WGS.shape[0])
print(overlap_imputed.shape[0])
overlap_id=pd.read_table(f"../WGS_IMP_Overlap/"+chrom+"/sites.txt",header=None).iloc[:, 0:4].astype(str).agg(':'.join, axis=1)

def process_row(i):
	a = overlap_imputed.iloc[i]
	b = overlap_WGS.iloc[i]
	return get_metrics(a, b,n_samp)

# Multiprocessing pool with tqdm progress bar
with mp.Pool(processes=mp.cpu_count()) as pool:
	results = list(tqdm(pool.imap(process_row, range(overlap_id.shape[0])), total=overlap_id.shape[0]))

# Combine the results into a DataFrame
out = pd.concat(results, axis=1).T.reset_index(drop=True)
out.columns = ['CHROM','POS','ID','REF','ALT','TP_HomR','TP_Het','TP_HomA','FN_HomR','FN_Het','FN_HomA','FP_HomR','FP_Het','FP_HomA','TN_HomR','TN_Het','TN_HomA']
print(out.head())
out.to_csv(f"{chrom}/UKB_7628_9768_ConfusionMetrics_chr{chrom}_SNPS_INDELs_Schunk{str(chunk_id)}.csv", index=False)
print("Done.\nTime elapsed for chr{chrom} chunk {str(chunk_id)}=",time.time()-start)


