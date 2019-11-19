# HF_GWAS
Heart Failure GWAS scripts and summary data

This repository contains example scripts used for the GWAS meta-analysis in the following paper:

Genome-wide association and multi-omic analyses reveal new mechanisms for Heart Failure
Marios Arvanitis, Yanxiao Zhang, Wei Wang, Adam Auton, 23AndMe Research Team, Ali R Keramati, Neil C. Chi, Bing Ren, Wendy S. Post, Alexis Battle
medRxiv 19006510; doi: https://doi.org/10.1101/19006510

The scripts require that the user has access to genotype files for a large population cohort in plink binary format (bed/bim/fam) 
in addition to phenotype files that determine at the least age, sex, race and the binary phenotype of interest for the genotyped individuals.

The user must also have access to a snpBatch file with SNP rsids and coordinates for the genotyped SNPs. 
That file is usually provided by the company doing the genotyping but if not available 
one can use the 1000Genomes snpBatch file instead which is available to download here: 
https://www.ncbi.nlm.nih.gov/projects/SNP/snp_viewBatch.cgi?sbid=1061891.

The scripts also require that the following executables are in your path:
1. plink (v1.9)
2. plink2 (v2a1)
3. bcftools
4. liftover (UCSC liftover tool) with the appropriate chain files
5. R (v3.4 or above)
6. python (v2.7) with pandas, subprocess and numpy installed
7. perl

The following R libraries are required:
1. argparser
2. data.table
3. SNPRelate
4. plyr
5. lattice
6. optparse
7. SAIGE

The pipeline also requires that the user downloads the following scripts from the supplemental data in the paper:
Anderson CA, Pettersson FH, Clarke GM, Cardon LR, Morris AP, Zondervan KT. 
Data quality control in genetic case-control association studies. 
Nat Protoc. 2010;5(9):1564–1573. doi:10.1038/nprot.2010.116
1. lmiss-hist.Rscript
2. run-diffmiss-qc.pl
3. run-IBD-QC.pl
*Note the first of those scripts should be edited to take input and output argparser arguments.

The GWAS steps consist of the following:

Step 1.
Run the Initial_genotype_qc.sh file which performs quality control filtering steps for the genotypes at the SNP and individual level.

Step 2.
Upload the output vcf files to the Michigan Imputation Server to get imputed genotypes

Step 3.
Download the imputed vcf files and run the script: download_and_analyze_imputed_data.sh to perform the GWAS in each cohort

Step 4.
Use METAL to meta-analyze the output summary statistics for each cohort.

The perl script for final QC checks before imputation was created by the McCarthy group and modified as appropriate:
https://www.well.ox.ac.uk/~wrayner/tools/#Checking

The R scripts to run SAIGE (step1_fitNULLGLMM.R and step2_SPAtests.R) were obtained as part of the SAIGE R package (https://github.com/weizhouUMICH/SAIGE)

