library(data.table)

data<-fread("unrelated_pfile.psam")

fwrite(data[,2], "indslist", quote=F, row.names=F, col.names=F, sep="\t", na="NA")

q()
n


