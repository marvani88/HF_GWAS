import argparse
import subprocess
import pandas as pd
import numpy as np

def check_bim_names(bim, snpBatch):
    """Read plink binary variant file names and assess format"""
    data=pd.read_csv(bim, sep='\t', header=None, dtype={'0':str, '1':str, '2':int, '3':int, '4':str, '5':str})
    data1=data.loc[data[1].str.contains('rs')]
    if len(data1)>0.75*len(data):
        print 'rsid'
        return 'rsid'
    else:
        data1=data.loc[data[1].str.contains('ss')]
        if len(data1)>0.75*len(data):
            print 'ss'
            return 'ss'
        else:
            data1=data[1][(len(data)/2)]
            data3=data[1][0]
            print data1
            print data3
            data2=pd.read_csv(snpBatch, skiprows=29, sep='\t', dtype={'#ss#':str, 'loc_snp_id':str, 'allele':str, 'samplesize':int, 'rs#':str, 'ss2rs_orien':str, 'chr':str, 'chr_pos':str, 'contig_acc':str, 'contig_pos':str, 'rs2genome_orien':str, 'assembly':str, 'weight':str})
            if data2['loc_snp_id'].str.contains(data1).any() and data2['loc_snp_id'].str.contains(data3).any():
                print 'locid'
                return 'locid'
            else:
                raise ValueError("Bim snp names not found in snpBatch")

def create_snp_map(scriptsdir, bim, snpBatch, indicator):
    """Change SNP coordinates to match those on snpBatch file if that is available and exclude the SNPs whose names are not found in the file"""
    indicator=check_bim_names(bim, snpBatch)
    if indicator=='rsid':
        subprocess.call(['Rscript', scriptsdir+'update_snp_map.R', snpBatch, '--bim', bim, '--rsid', 'rsid'])
        return True
    elif indicator=='locid':
        subprocess.call(['Rscript', scriptsdir+'update_snp_map.R', snpBatch, '--bim', bim, '--rsid', 'locid'])
    elif indicator=='ss':
        subprocess.call(['Rscript', scriptsdir+'update_snp_map.R', snpBatch, '--bim', bim, '--rsid', 'ss'])
    else:
        raise ValueError("Not an appropriate indicator value")

def liftover_snp_names_and_ids(scriptsdir,bim, snpBatch, genostart, hg, force):
    """Liftover coordinates to hg38 and remove unlifted SNPs:
    If >75% of the provided SNPBatch file SNP names matched the provided variant names we use that file for the liftover,
    otherwise liftover is performed based on provided SNP coordinates and the UCSC liftover tool"""
    file=genostart+'_merged'
    outrsid=genostart+'_rsid'
    outlift=genostart+'_old_hg38'
    outhg=genostart+'_hg38'
    print force
    if force=='yes':
        print "Forcing program to perform UCSC liftover"
        gencoords=open('genomicCoordinates.txt', 'w')
        GRC='/home-3/marvani1@jhu.edu/work/marios/UCSC_tools/chain_files/'+hg+'ToHg38.over.chain.gz'
        subprocess.call(["awk", '{print "chr"$1, $4-1, $4, $2}', bim], stdout=gencoords)
	coords=open('testcoords', 'w')
	subprocess.call(["awk", '{if ($1=="chr23") $1="chrX"}1', "genomicCoordinates.txt"], stdout=coords)
	subprocess.call(['mv testcoords genomicCoordinates.txt'], shell=True)
        subprocess.call(['liftOver', 'genomicCoordinates.txt', GRC, 'genomicCoordinates_hg38.txt', 'unlifted.bed'])
        subprocess.call(['Rscript', scriptsdir+'update_snp_map.R', snpBatch, '--bim', 'genomicCoordinates_hg38.txt', '--snppos', 'TRUE'])
        subprocess.call(['plink', '--bfile', file, '--update-map', 'snp.map', '4', '1', '--update-chr', 'snp.map', '3', '1', '--exclude', 'snps_to_remove', '--make-bed', '--out', outlift])
        subprocess.call(['plink', '--bfile', outlift, '--update-name', 'snp.map', '2', '1', '--make-bed', '--out', outhg])
    else:
        try:
            indicator=check_bim_names(bim, snpBatch)
            x=create_snp_map(bim, snpBatch, indicator)
            if x==True:
                subprocess.call(['plink', '--bfile', file, '--update-map', 'snp.map', '3', '1', '--update-chr', 'snp.map', '2', '1', '--exclude', 'snps_to_remove', '--make-bed', '--out', outhg])
            else:
                subprocess.call(['plink', '--bfile', file, '--update-name', 'snp.map', '2', '1', '--exclude', 'snps_to_remove', '--make-bed', '--out', outrsid])
                subprocess.call(['plink', '--bfile', outrsid, '--update-map', 'snp.map', '4', '2', '--update-chr', 'snp.map', '3', '2', '--make-bed', '--out', outhg])
        except ValueError:
            gencoords=open('genomicCoordinates.txt', 'w')
            GRC='/home-3/marvani1@jhu.edu/work/marios/UCSC_tools/chain_files/'+hg+'ToHg38.over.chain.gz'
            subprocess.call(["awk", '{print "chr"$1, $4-1, $4, $2}', bim], stdout=gencoords)
	    coords=open('testcoords', 'w')
            subprocess.call(["awk", '{if ($1=="chr23") $1="chrX"}1', "genomicCoordinates.txt"], stdout=coords)
            subprocess.call(['mv testcoords genomicCoordinates.txt'], shell=True)
            subprocess.call(['liftOver', 'genomicCoordinates.txt', GRC, 'genomicCoordinates_hg38.txt', 'unlifted.bed'])
            subprocess.call(['Rscript', scriptsdir+'update_snp_map.R', snpBatch, '--bim', 'genomicCoordinates_hg38.txt', '--snppos', 'TRUE'])
            subprocess.call(['plink', '--bfile', file, '--update-map', 'snp.map', '4', '1', '--update-chr', 'snp.map', '3', '1', '--exclude', 'snps_to_remove', '--make-bed', '--out', outlift])
            subprocess.call(['plink', '--bfile', outlift, '--update-name', 'snp.map', '2', '1', '--make-bed', '--out', outhg])

def clean_hg38_data(genostart):
    input=genostart+'_hg38'
    out=genostart+'_split_hg38'
    subprocess.call(['plink', '--bfile', input, '--split-x', 'b38','no-fail', '--make-bed', '--out', out])
    subprocess.call(['plink', '--bfile', out, '--make-bed', '--out', input])

def sex_check(scriptsdir, genostart, hg):
    """Check to see if self reported sex matches genotype sex and filter individuals in which it doesn't"""
    input=genostart+'_merged'
    if hg=='hg19':
        try:
            subprocess.check_call(['plink', '--bfile', input, '--split-x', 'b37', '--make-bed', '--out', 'splitx'])
            splitx='splitx'
        except subprocess.CalledProcessError:
            splitx=input
    elif hg=='hg18':
        try:
            subprocess.check_call(['plink', '--bfile', input, '--split-x', 'b36', '--make-bed', '--out', 'splitx'])
            splitx='splitx'
        except subprocess.CalledProcessError:
            splitx=input
    else:
        splitx=input
    subprocess.call(['plink', '--bfile', splitx, '--check-sex', '0.4', '0.8', '--out', 'sexstat'])
    subprocess.call(['Rscript', scriptsdir+'sex_check_plot.R', './sexstat.sexcheck', '--output', './sex_het_hist.pdf'])
    f=open('sex_check_fail.set', 'w')
    subprocess.call(['grep', 'PROBLEM', 'sexstat.sexcheck'], stdout=f)
    k=open('fail_sex_check.txt', 'w')
    subprocess.call(['awk', "{print $1, $2}", 'sex_check_fail.set'], stdout=k)

def het_miss(scriptsdir, genostart):
    """Filter individuals with extremes of heterozygosity and missingness"""
    input=genostart+'_hg38'
    subprocess.call(['plink', '--bfile', input, '--missing', '--out', 'miss'])
    subprocess.call(['plink', '--bfile', input, '--het', '--out', 'het'])
    subprocess.call(['Rscript', scriptsdir+'miss_het_plot.R', '--het', './het.het', '--imiss', './miss.imiss', '--plot', './gencall_het_plot.pdf', '--output', './fail_imiss_het.txt'])

def eliminate_dups(scriptsdir, genostart):
    """Filter duplicate samples"""
    subprocess.call(['cp', 'miss.imiss', 'dups.imiss'])
    subprocess.call(['cp', 'dups', 'dups.genome'])
    subprocess.call(['perl', scriptsdir+'run-IBD-QC.pl', 'dups'])

def clean_inds(genostart):
    """Exclude all individuals that failed the filters above.
    In addition, we exclude individuals whose self-reported ancestry does not match the 1000G europeans ancestry based on the first two principal components in a joint analysis
    (inds_to_remove.txt)-Script for that exclusion is not provided as it requires manually visualizing the joint PCA plots and choosing appropriate exclusion thresholds"""
    input=genostart+'_hg38'
    output='clean_inds_'+ genostart
    out=open('fail-qc-inds.txt', 'w')
    subprocess.call(['cat fail_imiss_het.txt fail_sex_check.txt fail-IBD-QC.txt inds_to_remove.txt| sort -k1 | uniq'], stdout=out, shell=True)
    subprocess.call(['plink', '--bfile', input, '--remove', 'fail-qc-inds.txt', '--make-bed', '--out', output])

def snp_call_rate_plot(scriptsdir, genostart):
    """Make a histogram plot with the SNP call rate"""
    input='clean_inds_'+ genostart
    subprocess.call(['plink', '--bfile', input, '--missing', '--out', 'lmiss'])
    subprocess.call(['Rscript', scriptsdir+'lmiss-hist.Rscript', '--lmiss', './lmiss.lmiss', '--output', './lmiss_hist.pdf'])

def snp_diffmiss(scriptsdir,genostart, pheno):
    """Assess differential missingness in HF cases and controls"""
    input='clean_inds_'+ genostart
    subprocess.call(['plink', '--bfile', input, '--pheno', pheno, '--test-missing', '--out', input, '--1'])
    subprocess.call(['perl', scriptsdir+'run-diffmiss-qc.pl', input])

def hardy(genostart, pheno, race, rno):
    """Exclude SNPs that have an extreme HWE deviaiton p-value (<1e-4) in any ancestry"""
    input='clean_inds_'+ genostart
    subprocess.call(['plink', '--bfile', input, '--pheno', pheno, '--1', '--filter-controls', '--make-bed', '--out', 'controls_only'])
    for num in range(rno):
        raceno=str(num+1)
        failfile='fail_hardy'+raceno+'.txt'
        subprocess.call(['plink', '--bfile', 'controls_only', '--geno', '0.05', '--pheno', race, '--mpheno', raceno, '--1', '--hardy', '--out', 'hardy'])
        newhardy=open('newhardy', 'w')
        subprocess.call(['grep', 'UNAFF', 'hardy.hwe'], stdout=newhardy)
        outfile=open(failfile, 'w')
        subprocess.call(['awk', "{if ($9<0.0001) print $2}", 'newhardy'], stdout=outfile)
    subprocess.call(['cat fail-diffmiss-qc.txt genotype_diffs fail_hardy* | sort -k1 | uniq > fail_qc_snps.txt'], shell=True)

def clean_snps(genostart, pheno):
    """Exclude SNPs that failed HWE or differential missingness filters or SNPs with MAF<0.01 or call rate <95%"""
    input='clean_inds_'+ genostart
    output='clean_' + genostart + '_data'
    subprocess.call(['plink', '--bfile', input, '--pheno', pheno, '--1', '--exclude', 'fail_qc_snps.txt', '--maf', '0.01', '--geno', '0.05', '--make-bed', '--out', output])

def lift_to_hg19(scriptsdir,genostart):
    """Liftover coordinates to hg19"""
    bim='final_' + genostart + '_data.bim'
    outlift='final_'+genostart+'_hg19'
    file='final_' + genostart + '_data'
    gencoords=open('genomicCoordinates_old38.txt', 'w')
    GRC='/home-3/marvani1@jhu.edu/work/marios/UCSC_tools/chain_files/hg38ToHg19.over.chain.gz'
    subprocess.call(["awk", '{print "chr"$1, $4-1, $4, $2}', bim], stdout=gencoords)
    coords=open('testcoords', 'w')
    subprocess.call(["awk", '{if ($1=="chr23") $1="chrX"}1', "genomicCoordinates_old38.txt"], stdout=coords)
    subprocess.call(['mv testcoords genomicCoordinates_old38.txt'], shell=True)
    subprocess.call(['liftOver', 'genomicCoordinates_old38.txt', GRC, 'genomicCoordinates_hg19.txt', 'unlifted_hg19.bed'])
    subprocess.call(['Rscript', scriptsdir+'snp_map_hg19.R'])
    subprocess.call(['plink', '--bfile', file, '--update-map', 'genomicCoordinates_hg19.snpmap', '3', '1', '--update-chr', 'genomicCoordinates_hg19.snpmap', '2', '1', '--allow-extra-chr', '--exclude', 'unlifted_hg19.snps', '--make-bed', '--out', 'test_hg19'])
    subprocess.call(['plink', '--bfile', 'test_hg19', '--allow-extra-chr', '--chr', '1-23', 'X', '--make-bed', '--out', outlift])
    subprocess.call(['rm test_hg19*'], shell=True)

def genotype_qc_basic(scriptsdir, genostart, pheno, snpBatch, race, rno, hg='hg19', force='no'):
    """Function incorporating all basic filtering steps described above"""
    bim=genostart+'_merged.bim'
    endfile='clean_' + genostart + '_data'
    finfile='final_' + genostart + '_data'
    liftover_snp_names_and_ids(scriptsdir,bim, snpBatch, genostart, hg, force)
    clean_hg38_data(genostart)
    sex_check(scriptsdir,genostart, hg)
    het_miss(scriptsdir,genostart)
    eliminate_dups(scriptsdir,genostart)
    clean_inds(genostart)
    snp_call_rate_plot(scriptsdir,genostart)
    snp_diffmiss(scriptsdir,genostart, pheno)
    hardy(genostart, pheno, race, rno)
    clean_snps(genostart, pheno)
    subprocess.call(['plink', '--bfile', endfile, '--merge-x', 'no-fail', '--make-bed', '--out', finfile])
    lift_to_hg19(scriptsdir,genostart)

##==========================Main===============================##

parser = argparse.ArgumentParser(description='Genotype_qc.')
parser.add_argument('--genostart', type=str,
                    help='string that contains the name of the merged genotypes')
parser.add_argument('--pheno', type=str,
                    help='phenotype fam file for plink with its corresponding path')
parser.add_argument('--snpBatch', default=False,
                    help='snpBatch file to use for rsid and chr:pos liftover')
parser.add_argument('--race', type=str,
                    help='race fam file for plink with its corresponding path')
parser.add_argument('--rno', type=int,
                    help='number of races in the study')
parser.add_argument('--force', default='no',
                    help='Do you want to force the program to use UCSC liftOver instead of matching based on snpBatch?')
parser.add_argument('--hg', default='hg19',
                    help='what is the base genome build')
parser.add_argument('--scriptsdir', default='./',
                    help='what is the directory of the scripts')

args = parser.parse_args()

subprocess.call(['rm fail* sexstat.sexcheck'], shell=True)
genotype_qc_basic(args.scriptsdir, args.genostart, args.pheno, args.snpBatch, args.race, args.rno, args.hg, args.force)
