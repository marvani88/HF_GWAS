library(data.table)
library(argparser)

p<-arg_parser("Update PCA columns on phenotype file")
p<-add_argument(p, "--pheno", help="Input a phenotype file", default="pheno.SAIGE.txt")
p<-add_argument(p, "--eigen", help="Input a file with new eigenvectors", default="")
p<-add_argument(p, "--out", help="Output phenotype file", default="")
argv<-parse_args(p)

data<-fread(argv$pheno, header=T)
ids<-data$IID

if("PC1" %in% colnames(data)){
x<-which(colnames(data)=="PC1")
data<-data[,1:(x-1)]
}

eigen<-fread(argv$eigen, header=T)

if("FID" %in% colnames(data)){
colnames(eigen)[1]<-"FID"
pheno<-merge(data, eigen, by="FID")
}else{
colnames(eigen)[1]<-"IID"
pheno<-merge(data, eigen, by="IID")
}

#if(dim(pheno)[1]!=dim(data)[1]){stop("Dimensions differ between data and pheno")}

pheno<-pheno[order(match(pheno$IID, ids)),]
fwrite(pheno, argv$out, quote=F, row.names=F, sep="\t", na="NA")

q()
n

