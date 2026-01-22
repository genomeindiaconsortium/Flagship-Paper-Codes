#!/bin/bash
# Separate GI samples into populations
mkdir -p populations
awk -F '\t' -v OFS='\t' '{print $1, $2 > "populations/"$4".txt"}' pop_info.txt

# Calculating frequency of variants per population
header=false
mkdir -p pop_freq_alt
for file in populations/*.txt;
do
	population=$(basename "$file" .txt)
	plink --bfile pharma_variants --keep "$file" --a2-allele ref_allele_fix.txt --freq --out "pop_freq_alt/$population"
done

# Fisher's exact test comparing populations
for file in populations/*.txt;
do
	population_id=$(basename "$file" .txt)
	mkdir $population_id
	plink --bfile pharma_variants --keep "$file" --a2-allele ref_allele_fix.txt --make-bed --out "$population_id/$population_id"
    plink --bfile pharma_variants --remove "$file" --a2-allele ref_allele_fix.txt --make-bed --out "$population_id/${population_id}_removed"
	awk '{$6=2;print}' "$population_id/${population_id}.fam" > "$population_id/${population_id}_cases.fam"
	mv "$population_id/${population_id}_cases.fam" "$population_id/${population_id}.fam"
	awk '{$6=1;print}' "$population_id/${population_id}_removed.fam" > "$population_id/${population_id}_controls.fam"
	mv "$population_id/${population_id}_controls.fam" "$population_id/${population_id}_removed.fam"
	plink --bfile "$population_id/${population_id}" --bmerge "$population_id/${population_id}_removed.bed" "$population_id/${population_id}_removed.bim" "$population_id/${population_id}_removed.fam" --assoc --allow-no-sex --a2-allele ref_allele_fix.txt --out "$population_id/${population_id}_fisher"
	if [ "$header" = false ]; then
		awk -v pop="$population_id" 'NR==1{if(FNR==1) print $0, "Populations"} NR>1 {print $0, pop}' "$population_id/${population_id}_fisher.assoc" > combined_assoc.assoc
	header=true
	else
		awk -v pop="$population_id" 'NR>1 {print $0, pop}' "$population_id/${population_id}_fisher.assoc" >> combined_assoc.assoc
	fi
done
