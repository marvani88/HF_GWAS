library(SNPRelate)
library(plyr)
library(data.table)
library(argparser)

p<-arg_parser("Remove related using SNPRelate")
p<-add_argument(p, "--vcf", help="Input a genotype VCF file", default="")
p<-add_argument(p, "--bed", help="Input a prefix for bed/bim/fam file", default="")
p<-add_argument(p, "--gds", help="Input a genotype GDS file", default="runPCA.gds")
p<-add_argument(p, "--prune", help="Do you want to do ld pruning?", default="no")
p<-add_argument(p, "--out", help="Output eigenvectors file", default=tempfile())
p<-add_argument(p, "--savegds", help="Do you want to save the gds file?", default=TRUE)
p<-add_argument(p, "--ibdthres", help="ibd threshold", type="float", default=0.125)
p<-add_argument(p, "--threads", help="Number of threads", type="int", default=1)
argv<-parse_args(p)

if(argv$vcf!=""){
  snpgdsVCF2GDS(argv$geno, argv$gds)}

if(argv$bed!=""){
  bed=paste(argv$bed, ".bed", sep="")
  bim=paste(argv$bed, ".bim", sep="")
  fam=paste(argv$bed, ".fam", sep="")
  snpgdsBED2GDS(bed, fam, bim, argv$gds)
}

genofile <- snpgdsOpen(argv$gds)

if(argv$prune!="no"){
  snpset<-snpgdsLDpruning(genofile, method="corr", slide.max.bp = 500000, ld.threshold = 0.1)
  snp.id <- unlist(snpset)}else{snp.id=NULL}

samples<-read.gdsn(index.gdsn(genofile, "sample.id"))

ibd.robust <- snpgdsIBDKING(genofile, snp.id=snp.id, num.thread=as.numeric(argv$threads))
sel<-snpgdsIBDSelection(ibd.robust, kinship=argv$ibdthres)
related<-NULL
while(nrow(sel)>0){
  sample.counts<-arrange(count(c(sel$ID1, sel$ID2)), -freq)
  rm.sample<-sample.counts[1, 'x']
  sel<-sel[sel$ID1!=rm.sample & sel$ID2!=rm.sample,]
  related<-c(as.character(rm.sample), related)
}
samples<-subset(samples, !(samples %in% related))

samples<-as.data.table(samples)

samples$V2<-samples$samples

fwrite(samples, argv$out, quote=F, row.names=F, col.names=F, sep="\t")

snpgdsClose(genofile)

if(argv$savegds!=TRUE){unlink("runPCA.gds", force=TRUE)}
