library(SNPRelate)
library(plyr)
library(data.table)
library(argparser)

p<-arg_parser("Perform PCA using SNPRelate")
p<-add_argument(p, "--vcf", help="Input a genotype VCF file", default="")
p<-add_argument(p, "--bed", help="Input a prefix for bed/bim/fam file", default="")
p<-add_argument(p, "--gds", help="Input a genotype GDS file", default="runPCA.gds")
p<-add_argument(p, "--prune", help="Do you want to do ld pruning?", default="no")
p<-add_argument(p, "--ibd", help="Do you want to runPCA only on unrelated samples?", default="no")
p<-add_argument(p, "--pdfpca", help="Output pdf pca file", default=tempfile())
p<-add_argument(p, "--pdfcr", help="Output pdf snp/pca correlation file", default="")
p<-add_argument(p, "--out", help="Output eigenvectors file", default=tempfile())
p<-add_argument(p, "--colorrace", help="Do you want to color plot by race?", default="no")
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

if(argv$ibd!="no"){
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
}

pca<-snpgdsPCA(genofile, sample.id=samples, snp.id=snp.id, num.thread=as.numeric(argv$threads))

if(argv$ibd!="no"){
  load<-snpgdsPCASNPLoading(pca, genofile, num.thread=as.numeric(argv$threads))
  finpca<-snpgdsPCASampLoading(load, genofile, num.thread=as.numeric(argv$threads))}else{finpca<-pca}

var<-data.frame(EIG=1:length(pca$varprop), Proportion=round(pca$varprop*100, 2))

if(argv$colorrace!="no"){
  pop_code <- read.gdsn(index.gdsn(genofile, "sample.annot"))
  sample<-read.gdsn(index.gdsn(genofile, "sample.id"))
  pop = factor(pop_code$phenotype)[match(finpca$sample.id, sample)]
}else{pop=1}

tab <- data.table(sample.id = finpca$sample.id, finpca$eigenvect[,1:10], pop=pop, stringsAsFactors=F)
c<-NULL
for(i in seq(10)){c<-c(c, paste("PC", i, sep=""))}
colnames(tab)<-c("samples",c, "race")
cab<-tab[,c(1:11)]
fwrite(cab, argv$out, quote=F, row.names=F, sep="\t")

chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))

png(argv$pdfpca)
par(mfrow=c(2,3))
plot(var$EIG, var$Proportion, xlab="Principal Component",ylab="Variance explained (%)", xlim=c(0,10))
plot(tab$PC1, tab$PC2,col=as.integer(tab$race))
legend("bottomright", legend=levels(tab$race), pch="o", col=1:nlevels(tab$race))
plot(tab$PC3, tab$PC4, col=as.integer(tab$race))
legend("bottomright", legend=levels(tab$race), pch="o", col=1:nlevels(tab$race))
plot(tab$PC5, tab$PC6, col=as.integer(tab$race))
legend("bottomright", legend=levels(tab$race), pch="o", col=1:nlevels(tab$race))
plot(tab$PC7, tab$PC8, col=as.integer(tab$race))
legend("bottomright", legend=levels(tab$race), pch="o", col=1:nlevels(tab$race))
plot(tab$PC9, tab$PC10, col=as.integer(tab$race))
legend("bottomright", legend=levels(tab$race), pch="o", col=1:nlevels(tab$race))
dev.off()

if(argv$pdfcr!=""){
  cr <- snpgdsPCACorr(pca, genofile, eig.which=1:10, snp.id=snp.id, num.thread=as.numeric(argv$threads))
  png(argv$pdfcr)
  par(mfrow=c(2,5))
  for(i in seq(10)){
    plot(abs(cr$snpcorr[i,]), xlab="SNP Index", ylab=paste("PC ", i, sep=""), col=chr)
    }
  dev.off()
}

snpgdsClose(genofile)

if(argv$savegds!=TRUE){unlink("runPCA.gds", force=TRUE)}
