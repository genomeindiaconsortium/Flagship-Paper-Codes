import time
import os
from mpi4py import MPI
import re
import gzip


# Define the list of chromosomes using a more readable format
chromosomes = [f'{i}' for i in range(1, 23)]
chrom_num=[f'{i}' for i in range(1, 23)]

base_name="all_sorted_withID_MIE0.05Rem.vcf.gz"
mapdir='../Required/chr-maps'
fam_file='9768_Fam_file.txt'
thread_count = str(os.environ.get('SLURM_CPUS_PER_TASK', 1))

	
def common_phasing(num):
	chrom_dir = os.path.join(os.getcwd(), chromosomes[num])
	if not os.path.exists(chrom_dir):
                os.makedirs(chrom_dir)
	main_file=f'Input/9768_chr{chromosomes[num]}{base_name}'
	os.system(f"phase_common --input {main_file}  --filter-maf 0.001 --region chr{chrom_num[num]}  --map {mapdir}/chr{chromosomes[num]}.b38.gmap.gz --output {chrom_dir}/chr{chromosomes[num]}_{base_name}_common_phased.bcf  --thread {thread_count} --pedigree {fam_file}")
	os.system(f"time bcftools convert {chrom_dir}/chr{chromosomes[num]}_{base_name}_common_phased.bcf -Oz -o {chrom_dir}/chr{chromosomes[num]}_{base_name}_common_phased.vcf.gz --threads {thread_count}")
	#os.system(f"bgzip -@ {thread_count} {chrom_dir}/chr{chromosomes[num]}_{base_name}_common_phased.vcf")
	os.system(f"bcftools index -t --threads {thread_count} {chrom_dir}/chr{chromosomes[num]}_{base_name}_common_phased.vcf.gz")


def rare_phasing(num):
	chrom_dir = os.path.join(os.getcwd(), chromosomes[num])
	if not os.path.exists(chrom_dir):
		os.makedirs(chrom_dir)
	main_file=f'Input/9768_chr{chromosomes[num]}{base_name}'
	os.system(f"GLIMPSE2_chunk --input {main_file} --region chr{chrom_num[num]}  --map {mapdir}/chr{chromosomes[num]}.b38.gmap.gz --sequential --output {chrom_dir}/chunks_chr{chromosomes[num]}.txt")
	os.system(f"mkdir {chrom_dir}/chunks")
	with open(f'{chrom_dir}/chunks_chr{chromosomes[num]}.txt', 'r') as file:
		for line in file: 
			os.system(f'phase_rare --input {main_file} --scaffold {chrom_dir}/chr{chromosomes[num]}_{base_name}_common_phased.vcf.gz --map {mapdir}/chr{chromosomes[num]}.b38.gmap.gz --input-region {line.split()[3]} --scaffold-region {line.split()[2]} --output {chrom_dir}/chunks/chr{chromosomes[num]}_{base_name}.target.phased.chunk{line.split()[0]}.vcf --thread {thread_count} --pedigree {fam_file}')
	os.system(f"ls -1v  {chrom_dir}/chunks/chr{chromosomes[num]}_{base_name}.target.phased.chunk*.vcf >  {chrom_dir}/chunks/chr{chromosomes[num]}_chunks_list.txt")
	os.system(f"bcftools concat -Oz --threads {thread_count} -o {chrom_dir}/chr{chromosomes[num]}_{base_name}.target.phased.vcf.gz -f {chrom_dir}/chunks/chr{chromosomes[num]}_chunks_list.txt")
	os.system(f"bcftools index -t --threads {thread_count}  {chrom_dir}/chr{chromosomes[num]}_{base_name}.target.phased.vcf.gz")
 
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank < len(chromosomes):
	common_phasing(rank)
	rare_phasing(rank)