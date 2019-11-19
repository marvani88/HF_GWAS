###Takes as input a 

library(argparser)
library(data.table)

p<-arg_parser("Create snp map and rsid")
p<-add_argument(p, "snplist", help="Input a snp Batch file")
p<-add_argument(p, "--bim", help="Input a bim file")
p<-add_argument(p, "--rsid", help="Do you need to update snp names", default="rsid")
p<-add_argument(p, "--snppos", help="Do you want to use chr:pos to update snp names (to be used only if you have an unlifted.bed file in the folder)", default=FALSE)
p<-add_argument(p, "--snpmap", help="Output snp.map file", default=tempfile())
p<-add_argument(p, "--snpremove", help="Output snps to remove file", default=tempfile())
argv<-parse_args(p)

x<-argv$snplist
snplist<-fread(x, header=T, stringsAsFactors=F, colClasses=list(character=8), skip=10)
bim<-fread(argv$bim, header=F, stringsAsFactors=F)

if(argv$snppos){
  snplist<-snplist[,c(5,7,8)]
  snplist<-unique(snplist)
  snplist$chrpos<-paste(snplist$chr, snplist$chr_pos, sep=":")
  gencoords<-bim
  colnames(gencoords)<-c("chr", "junk", "chr_pos", "name")
  gencoords$chr<-gsub("chr", "", gencoords$chr)
  gencoords$chrpos<-paste(gencoords$chr, gencoords$chr_pos, sep=":")
  unlifted<-read.table("unlifted.bed", comment.char="#", header=F, stringsAsFactors=F)
  snplist<-subset(snplist, snplist$chrpos %in% gencoords$chrpos)
  rem<-gencoords$name[which(!(gencoords$chrpos %in% snplist$chrpos))]
  snplist<-merge(snplist, gencoords, by="chrpos", all.x=T, all.y=T)
  to_remove<-snplist$name[which(snplist$chr_pos.x=="N.D."|snplist$chr.x=="N.D."|is.na(snplist$chr_pos.x))]
  snplist<-snplist[-which(snplist$chr_pos.x=="N.D."|snplist$chr.x=="N.D."|is.na(snplist$chr_pos.x)),]
  unlifted<-unlifted$V4
  snplist<-snplist[order(chr.x),]
  snplist<-snplist[!duplicated(snplist$name),]
  dups<-snplist$name[duplicated(snplist$rs)]
  snplist<-snplist[!duplicated(snplist$rs),]
  rem<-c(unlifted, rem, to_remove, dups)
  rem<-unique(rem)
  rem<-as.data.table(rem)
  snplist<-snplist[,c(8, 2, 3, 4)]
  print(paste("Number of snps in remove file: ", dim(rem)[1], sep=""))
  fwrite(rem, "snps_to_remove", quote=F, row.names=F, col.names=F, sep="\t")
  fwrite(snplist, "snp.map", quote=F, row.names=F, col.names=F, sep="\t")

}else if(argv$rsid=="locid"){
  snplist<-snplist[,c(2,5,7,8)]
  snplist$loc_snp_id<-gsub("AFFY_6_1M_", "", snplist$loc_snp_id)
  snplist<-unique(snplist)
  newlist<-bim$V2
  snplist<-subset(snplist, snplist$loc_snp_id %in% newlist)
  snplist<-snplist[order(chr),]
  snplist<-snplist[!duplicated(snplist$rs),]
  to_remove<-snplist$loc_snp_id[which(snplist$chr_pos=="N.D."|snplist$chr=="N.D.")]
  snplist<-snplist[-which(snplist$chr_pos=="N.D."|snplist$chr=="N.D."),]
  rem<-bim$V2[which(!(bim$V2 %in% snplist$loc_snp_id))]
  to_remove<-c(to_remove, rem)
  to_remove<-unique(to_remove)
  print(paste("Number of snps in remove file: ", length(to_remove), sep=""))
  to_remove<-as.data.table(to_remove)
  fwrite(to_remove, "snps_to_remove", quote=F, row.names=F, col.names=F, sep="\t")
  fwrite(snplist, "snp.map", quote=F, row.names=F, col.names=F, sep="\t")
  
  }else if(argv$rsid=="ss"){
    snplist<-snplist[,c(1,5,7,8)]
    snplist<-unique(snplist)
    newlist<-bim$V2
    colnames(snplist)[1]<-"loc_snp_id"
    snplist<-subset(snplist, snplist$loc_snp_id %in% newlist)
    snplist<-snplist[order(chr),]
    snplist<-snplist[!duplicated(snplist$rs),]
    to_remove<-snplist$loc_snp_id[which(snplist$chr_pos=="N.D."|snplist$chr=="N.D.")]
    snplist<-snplist[-which(snplist$chr_pos=="N.D."|snplist$chr=="N.D."),]
    rem<-bim$V2[which(!(bim$V2 %in% snplist$loc_snp_id))]
    to_remove<-c(to_remove, rem)
    to_remove<-unique(to_remove)
    print(paste("Number of snps in remove file: ", length(to_remove), sep=""))
    to_remove<-as.data.table(to_remove)
    fwrite(to_remove, "snps_to_remove", quote=F, row.names=F, col.names=F, sep="\t")
    fwrite(snplist, "snp.map", quote=F, row.names=F, col.names=F, sep="\t")
  
  }else{
  snplist<-snplist[,c(5,7,8)]
  snplist<-unique(snplist)
  newlist<-bim$V2
  snplist<-subset(snplist, snplist$rs %in% newlist)
  snplist<-snplist[order(chr),]
  snplist<-snplist[!duplicated(snplist$rs),]
  to_remove<-snplist$loc_snp_id[which(snplist$chr_pos=="N.D."|snplist$chr=="N.D.")]
  snplist<-snplist[-which(snplist$chr_pos=="N.D."|snplist$chr=="N.D."),]
  rem<-bim$V2[which(!(bim$V2 %in% snplist$rs))]
  to_remove<-c(to_remove, rem)
  to_remove<-unique(to_remove)
  print(paste("Number of snps in remove file: ", length(to_remove), sep=""))
  to_remove<-as.data.table(to_remove)
  fwrite(to_remove, "snps_to_remove", quote=F, row.names=F, col.names=F, sep="\t")
  fwrite(snplist, "snp.map", quote=F, row.names=F, col.names=F, sep="\t")}

q()
n
