library(data.table)
library(argparser)

p<-arg_parser("Transform SNP ID")
p<-add_argument(p, "--input", help="Input study")
p<-add_argument(p, "--out", help="Output study")
argv<-parse_args(p)

data<-fread(argv$input)
data[,"SNPID":=paste(CHR, ":", POS, "@", Allele1, "@", Allele2, sep="")]
fwrite(data, argv$out, quote=F, row.names=F, sep="\t")
q()
n


