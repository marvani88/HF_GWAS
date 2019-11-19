#.libPaths("~/work/marios/R_package/")
library(argparser)
library(data.table)

##Parse arguments
p<-arg_parser("Missingness and heterozygosity plot")
p<-add_argument(p, "--het", help="Input plink het file containing heterozygosity rate info")
p<-add_argument(p, "--imiss", help="Input plink het file containing individual missingness info")
p<-add_argument(p, "--plot", help="Output scatterplot file", default=tempfile())
p<-add_argument(p, "--output", help="Output failed individuals file", default=tempfile())
p<-add_argument(p, "--thres", help="Set the missingness threshold manually", default=-1)
argv<-parse_args(p)

##Merge plink data
het<-fread(argv$het, header=T, stringsAsFactors=F)
miss<-fread(argv$imiss, header=T, stringsAsFactors=F)
het$FHET<-(het$`N(NM)`-het$`O(HOM)`)/het$`N(NM)`
miss<-miss[,c(1,2,6)]
het<-het[,c(1,2,7)]
hetmiss<-merge(het, miss, by="IID")

##Calculate mean +/- 3 sd
h1=mean(hetmiss$FHET)+3*sd(hetmiss$FHET)
h2=mean(hetmiss$FHET)-3*sd(hetmiss$FHET)
v1=mean(hetmiss$F_MISS)+3*sd(hetmiss$F_MISS)
v2=mean(hetmiss$F_MISS)-3*sd(hetmiss$F_MISS)
if(argv$thres>=0){v=argv$thres}else{v=min(0.1,v1)}
if(v==0){v=0.003}

#if(length(subset(hetmiss$IID, !(hetmiss$F_MISS<v&hetmiss$FHET>h2&hetmiss$FHET<h1)))>0.9*length(hetmiss$IID)){v=max(v1, 0.03)}

print(paste("missingness cutoff:", v, sep=" "))
print(paste("Individuals removed for heterogeneity or missingness:", length(subset(hetmiss$IID, !(hetmiss$F_MISS<v&hetmiss$FHET>h2&hetmiss$FHET<h1)))), sep=" ")

##Make plot
pdf(argv$plot)
plot(hetmiss$F_MISS, hetmiss$FHET, type="p", col="blue", xlab="Proportion of missing genotypes", ylab="Heterozygosity rate", ylim=c(0.1, 0.7), xlim=c(0.001,0.5))
abline(h=h1, lty=2, col="red")
abline(h=h2, lty=2, col="red")
abline(v=v, lty=2, col="green")
dev.off()

hetmiss<-subset(hetmiss, !(hetmiss$F_MISS<v&hetmiss$FHET>h2&hetmiss$FHET<h1))
hetmiss<-hetmiss[,c(2,1)]

fwrite(hetmiss, argv$output, quote=F, sep="\t", row.names=F, col.names=F)
