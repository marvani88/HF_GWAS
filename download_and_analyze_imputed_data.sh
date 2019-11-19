#!/bin/bash -l
#SBATCH
#SBATCH --job-name=imputed_data
#SBATCH --partition=shared
#SBATCH --time=72:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=22
#SBATCH --mem=60GB
#SBATCH -o slurm.server.imputed.all.%j.out
#SBATCH -e slurm.server.imputed.all.%j.err
#SBATCH --account=abattle4

scriptsdir=$1
maindir=$2
phenotype_file=$3

cd ${maindir}

###Step 0 (not shown) >> Download and unzip imputed vcf files (one for each chromosome) by the Michigan Imputation Server

##Step 1 >> Filter R2>0.3 and merge files by chromosome

for chr in {1..22}
do
echo $chr
bcftools filter -i 'INFO/R2>0.3' chr${chr}.dose.vcf.gz -Oz -o TEMP_${chr}.vcf.gz &
done

wait

bcftools concat TEMP_{1..22}.vcf.gz -Oz -o study_filter_R2_0.3.vcf.gz --threads 21

##Step 2 >> filter MAF>0.01 and make plink2 binaries
bcftools filter -e 'INFO/MAF<0.01' study_filter_R2_0.3.vcf.gz -Oz -o study_filter_R2_0.3_MAF_0.01.vcf.gz --threads 21

tabix -p vcf study_filter_R2_0.3_MAF_0.01.vcf.gz && rm TEMP*

plink2 --vcf study_filter_R2_0.3_MAF_0.01.vcf.gz dosage=DS --out study --threads 22

##Step 3 >> filter MAF>0.01, missing call rate 0 (should not filter anything in imputed data), assess HWE p-value and filter HWE p<1e-4 to generate a pruned file for PCA and relatedness calculations 

plink2 --pfile study --pheno ${phenotype_file} --pheno-name chf --maf 0.01 --geno 0 --remove-if chf==1 --1 --hardy --out hardy --threads 22
awk '{if ($10<0.00001) print $2}' hardy.hardy > fail_hardy.txt

plink2 --pfile study --extract-if-info 'R2>0.9' --exclude fail_hardy.txt --geno 0 --maf 0.01 --indep-pairwise 500kb 1 0.2 --out pruned --threads 22

plink2 --pfile study --extract pruned.prune.in --make-bed --out hardcall_R2_0.9_ld_0.2_hwe --threads 22

mkdir unrelated

###Step 4 >> Remove all individuals that are related at an IBD threshold 0.125 and create plink binary file with unrelated individuals

Rscript ${scriptsdir}remove_related.R \
        --bed hardcall_R2_0.9_ld_0.2_hwe \
        --gds hardcall_ld_0.2_R2_0.9_hwe.gds \
        --threads 22 \
        --out ./unrelated/unrelated.samples 

plink2 --pfile study --extract-if-info 'R2>0.3' --keep ./unrelated/unrelated.samples --maf 0.01 --geno 0 --make-pgen --out ./unrelated/unrelated_pfile --threads 22

###Step 5 >> Filter HWE<1e-4 and generate high quality pruned SNPs for PCA calculation

plink2 --pfile unrelated/unrelated_pfile --pheno pheno_mendrand.txt --pheno-name CHF --keep-if CHF==control --1 --hardy --out unrelated/hardy_race1 --threads 22

awk '{if ($10<0.00001) print $2}' unrelated/hardy_race1.hardy > unrelated/fail_hardy.txt

cd unrelated

Rscript ${scriptsdir}filter_hardy.R --hardy hardy_race1.hardy

plink2 --pfile unrelated_pfile --pheno pheno_mendrand.txt --pheno-name CHF --require-pheno CHF --extract-if-info 'R2>0.9' --exclude fail_hardy.txt --indep-pairwise 500kb 1 0.2 --1 --out pruned --threads 22

plink2 --pfile unrelated_pfile --extract pruned.prune.in --pheno ${phenotype_file} --pheno-name CHF --require-pheno CHF --make-bed --out hardcall_R2_0.9_ld_0.2_hwe --1 --threads 22

###Step 6 >> Calculate principal components with a stricter relatedness threshold (0.08) and incorporate them as covariates in the phenotype file

Rscript ${scriptsdir}runPCA.R \
        --bed hardcall_R2_0.9_ld_0.2_hwe \
        --ibd yes \
	--ibdthres 0.08 \
        --pdfcr PC_chr_correlation_R2_0.9_ld_0.2_hwe.png \
        --gds hardcall_ld_0.2_R2_0.9_hwe.gds \
        --colorrace yes \
        --pdfpca PCA_R2_0.9_ld_0.2_hwe.png \
        --threads 22 \
        --out PCA_R2_0.9_ld_0.2_hwe_eigenvectors.tsv

Rscript ${scriptsdir}change_pheno_cols.R --pheno ${phenotype_file} --eigen PCA_R2_0.9_ld_0.2_hwe_eigenvectors.tsv --out ${phenotype_file}

###Step 7 >> Exclude SNPs that fail HWE filter and poorly imputed SNPs with R2<0.7 and reapply filter MAF>0.01 in unrelated individuals to generate the final VCFtools file for SAIGE

plink2 --pfile ${genotype_file} --maf 0.01 --geno 0.02 --extract-if-info 'R2>0.7' --exclude fail_hardy.txt --export vcf bgz vcf-dosage=DS-force --out SAIGE_unrelated --threads 22
tabix -p vcf SAIGE_unrelated.vcf.gz

Rscript ${scriptsdir}create_indslist.R

###Step 8 >> Run SAIGE with age, gender, 10 genotype PCs as covariates

mkdir SAIGE

Rscript ${scriptsdir}step1_fitNULLGLMM.R \
        --plinkFile=./hardcall_R2_0.9_ld_0.2_hwe \
        --phenoFile=${phenotype_file} \
        --phenoCol=CHF \
        --covarColList=age,gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
        --sampleIDColinphenoFile=IID \
        --traitType=binary \
        --outputPrefix=./SAIGE/SAIGE_CHF_all \
        --traceCVcutoff=0.0025 \
        --ratioCVcutoff=0.01 \
        --nThreads=22 \
        --LOCO=TRUE

for chr in {1..22}
do
echo chr${chr}
Rscript ${scriptsdir}step2_SPAtests.R \
        --vcfFile=./SAIGE_unrelated.vcf.gz \
        --vcfFileIndex=./SAIGE_unrelated.vcf.gz.tbi \
        --chrom=${chr} \
        --minMAC=3 \
        --minMAF=0.01 \
        --sampleFile=indslist \
        --GMMATmodelFile=./SAIGE/SAIGE_CHF_all.rda \
        --varianceRatioFile=./SAIGE/SAIGE_CHF_all.varianceRatio.txt \
        --SAIGEOutputFile=./SAIGE/SAIGE_CHF_all.${chr}.plainDosage.SAIGE.txt \
        --numLinesOutput=10000 \
        --IsOutputAFinCaseCtrl=TRUE \
        --LOCO=TRUE &
done

wait

head -1 ./SAIGE/SAIGE_CHF_all.1.plainDosage.SAIGE.txt > ./SAIGE/SAIGE_CHF_all.plainDosage.SAIGE.txt && tail -n +2 -q ./SAIGE/SAIGE_CHF_all.{1..22}.plainDosage.SAIGE.txt >> ./SAIGE/SAIGE_CHF_all.plainDosage.SAIGE.txt

rm ./SAIGE/SAIGE_CHF_all.{1..22}.plainDosage.SAIGE.txt

cd SAIGE

Rscript ${scriptsdir}SAIGE_plot.R --dosage SAIGE_CHF_all.plainDosage.SAIGE.txt --qqplot qqplot.SAIGE.CHF.png && rm ../SAIGE_unrelated.vcf*
