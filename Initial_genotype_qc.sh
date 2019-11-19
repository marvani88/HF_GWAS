#!/bin/bash -l
#SBATCH
#SBATCH --job-name=Substudy_plink_binaries_genotype_qc
#SBATCH --partition=lrgmem
#SBATCH --time=10:0:0
#SBATCH --nodes=1
#SBATCH --mem=100GB
#SBATCH -o slurm.Substudy_plink_binaries.genoqc.%j.out
#SBATCH -e slurm.Substudy_plink_binaries.genoqc.%j.err
#SBATCH --account=abattle4

scriptsdir=$1
maindir=$2
snpBatchfile=$3
genotypes=$4
phenotypes=$5
race_file=$6
genome_build=$7

cd maindir

python ${scriptsdir}Update_removed_inds_list.py --genostart ${genotypes} --snpBatch ${snpBatchfile} --pheno ${phenotypes} --race ${race_file} --rno 1 --hg ${genome_build}

mkdir perlfiles_for_server

for FILE in final_${genotypes}_hg19*
do
mv ${FILE} ./perlfiles_for_server/${FILE}
done

cd ./perlfiles_for_server

plink --bfile final_${genotypes}_hg19 --freq --out final_${genotypes}_hg19 --nonfounders

perl HRC-1000G-check-bim-NoReadKey.pl -b final_${genotypes}_hg19.bim -f final_${genotypes}_hg19.frq -r ../annotations/1000GP_Phase3_combined.legend -g -p EUR

bash Run-plink.sh

mkdir ../preimpute

cd ../preimpute

for chr in {1..23}
do
echo chr${chr}
plink2 --bfile ../perlfiles_server/final_${genotypes}_hg19-updated --chr ${chr} --ref-allele force ../perlfiles_server/Force-Allele1-final_${genotypes}_hg19-1000G.txt --recode vcf-4.2 bgz --out ${genotypes}_chr${chr} --real-ref-alleles
done



