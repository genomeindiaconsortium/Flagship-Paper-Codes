import subprocess
import gzip
import os
from mpi4py import MPI

performing for Chr1-22
x=[f'{i}' for i in range(1,23)]
y=[str(i) for i in range(1,23)]



base_dir='/gpfs/data/user/shreya/Crucial/GI_9768' 
#ref_dir=f"{base_dir}/3.Phasing/"   ## Where reference panel has been created
ref_dir=f"{base_dir}/MAC3_Panel/"   # Where reference panel has been created
map_dir = f"{base_dir}/Required/chr-maps"            ##
sanscog_dir=f"{base_dir}/UKB_Indian/UKB_GI/INPUT/"
imp_dir=f"{base_dir}/UKB_Indian/UKB_GI/OUTPUT/"


thread_count = str(os.environ.get('SLURM_CPUS_PER_TASK', 1))
def imputation(num):
	os.system(f"mkdir -p {imp_dir}{x[num]}/out")
	os.system(f"mkdir {imp_dir}{x[num]}/log")
	
	os.system(f"imp5Chunker_1.1.5_static --h {ref_dir}chr{x[num]}_all_sorted_withID_MIE0.05Rem.vcf.gz.phased.annot_MAC3.vcf.gz --g {sanscog_dir}/chr{x[num]}_ukb22418_IND_PAK_BAN_ID.hg38_chr.renamed.QCd.MAC2Final.target.phased.paneloverlap.vcf.gz  --window-size 5000000  --r chr{y[num]}  --o {imp_dir}{x[num]}/{x[num]}-coordinates.txt")
	buffer_region=[]
	imp_region=[]
	with open(f"{imp_dir}{x[num]}/{x[num]}-coordinates.txt",'r') as chunks:
		for line in chunks:
			cols=line.strip().split('\t')
			#Saving 3rd column Buffer Region
			buffer_region.append(cols[2])
			#Saving 4th column Imputation Region
			imp_region.append(cols[3])
	print("Buffer Region\n",buffer_region,"\n Imputation Regions\n ",imp_region)   
	## Imputation process split across multiple chunks according to the Chromosomal Map
	c=0
	while (c<len(buffer_region)):
		z=str(c+1)
		os.system(f"time impute5_1.1.5_static --m {map_dir}/chr{x[num]}.b38.gmap.gz --h {ref_dir}/chr{x[num]}_all_sorted_withID_MIE0.05Rem.vcf.gz.phased.annot_MAC3.vcf.gz --g {sanscog_dir}/chr{x[num]}_ukb22418_IND_PAK_BAN_ID.hg38_chr.renamed.QCd.MAC2Final.target.phased.paneloverlap.vcf.gz   --r {imp_region[c]}  --o {imp_dir}{x[num]}/out/{x[num]}-chunk-{z}-imputed_MAC3.vcf.gz --l {imp_dir}{x[num]}/log/{x[num]}-chunk-{z}.log --threads {thread_count} --out-gp-field --buffer-region {buffer_region[c]}")
		c+=1
	## Merging all the chunks in order
	os.system(f"ls -1v {imp_dir}{x[num]}/out/{x[num]}-chunk-*-imputed_MAC3.vcf.gz > {imp_dir}{x[num]}/{x[num]}-chunks_list.txt")
	
	## Concatenation
	os.system(f"bcftools concat {imp_dir}{x[num]}/out/{x[num]}-chunk-*-imputed_MAC3.vcf.gz -Oz -o {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-mac3.vcf.gz --threads {thread_count}")
	## Sorting
	os.system(f"bcftools sort {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-mac3.vcf.gz -Oz -o {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-sorted-mac3.vcf.gz  --temp-dir {imp_dir}{x[num]}/TEMP")

	os.system(f"bcftools index -t {x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-sorted-mac3.vcf.gz --threads {thread_count}")
	## Adding MAF in INFO field
	os.system(f"bcftools +fill-tags --threads {thread_count} {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-sorted-mac3.vcf.gz  -Oz -o {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-unfiltered-sorted-mac3-annot.vcf.gz -- -t MAF")
	os.system(f"bcftools index -t  {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-unfiltered-sorted-mac3-annot.vcf.gz --threads {thread_count} -f")
	## Deleting the intermediate chunk files
	#os.system(f"rm -rf {x[num]}/out/")
	## Extracting the Variant ID, Allele Frequency and R.Sq. Values for Metric File creation
	os.system(f"bcftools query -f '%ID\t%INFO/AF\t%INFO/INFO\n' {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-unfiltered-sorted-mac3-annot.vcf.gz > {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-GI-unfiltered-Metrics.txt")
	os.system(f"bgzip -@ {thread_count} -f {imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-GI-unfiltered-Metrics.txt")



def postimp(num):
	main_file=f'{imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-imputed-GI-unfiltered-sorted-mac3-annot.vcf.gz'
	##Getting sample count for MAC Calculation
	command=f"bcftools query -l {main_file} |wc -l"
	sample_num= subprocess.run([command],text=True,stdout=subprocess.PIPE,shell=True,stderr=subprocess.PIPE)
	sample_count=float(sample_num.stdout)
#	print(sample_count)
	## Calculating MAC and MAF using python built in commands
	with gzip.open(f"{imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-GI-unfiltered-Metrics.txt.gz",'rt') as ip, gzip.open(f"{imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-GI-unfiltered-MAF-MAC-Metrics.txt.gz",'wt') as op_mac:
		for line in ip:
			cols=line.strip().split('\t')
			if float(cols[1]) > 0.5:
				maf_new=float(1-float(cols[1]))
				cols[1]=str(maf_new)
			mac_new = round((float(cols[1]) * (2*sample_count)))
			output = f"{cols[0]}\t{cols[1]}\t{mac_new:.0f}\t{cols[2]}\n"
			op_mac.write(output)
	metrics_file= gzip.open(f'{imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-GI-unfiltered-MAF-MAC-Metrics.txt.gz','rt')
	op6maf2ids= gzip.open(f'{imp_dir}{x[num]}/UKB_IND_PAK_BAN-{x[num]}-GI_info0.6-flip-0.0002_IDs.txt.gz','wt')
	next(metrics_file)
	for line in metrics_file:
		cols=line.strip().split('\t')
		id=cols[0].strip().split(':')
		ref=id[2]
		alt=id[3]
		maf=float(cols[1])
		info=float(cols[3])
		##The following is the criteria for Allele Flip(ATGC)
		flip_check = (0.4 <= maf <= 0.6) and ((ref in ['A','T'] and alt in ['A','T']) or (ref in ['G','C'] and alt in ['G','C']))
		if (info >= 0.6 and maf>=0.0002 and not flip_check):
			op6maf2ids.write(str(id[0])+'\t'+str(id[1])+'\n')
	

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank < (len(x)):
	imputation(rank)
	postimp(rank)
