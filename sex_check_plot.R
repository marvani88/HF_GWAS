#.libPaths("~/work/marios/R_package/")
library(argparser)
library(data.table)

##Parse arguments
p<-arg_parser("Sex heterozygosity histogram")
p<-add_argument(p, "input", help="Input plink sexcheck file containing sex heterozygosity info")
p<-add_argument(p, "--output", help="Output histogram file", default=tempfile())
argv<-parse_args(p)

##load data
mydata<-fread(argv$input, header=T, stringsAsFactors=F)

##Make plot
pdf(argv$output, width=6, height=6)
par(mfrow=c(2,1))
hist(subset(mydata$F, mydata$PEDSEX==2), freq=F, breaks=100, xlab="X chromosome homozygocity", ylab="Frequency", main="All female samples")
if(1 %in% mydata$PEDSEX){hist(subset(mydata$F, mydata$PEDSEX==1), freq=F, breaks=100, xlab="X chromosome homozygocity", ylab="Frequency", main="All male samples")}
dev.off()
