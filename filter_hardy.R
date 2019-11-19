library(data.table)
library(argparser)

p<-arg_parser("Filter based on hardy-weinberg threshold")
p<-add_argument(p, "--hardy", help="Input hardy weinberg data")
p<-add_argument(p, "--out", help="output rsids to be filtered", default="filter_hardy.txt")
argv<-parse_args(p)

data<-fread(argv$hardy)
data<-subset(data, data$P<1e-4)
data<-data[,2]
fwrite(data, "fail_hardy.txt", quote=F, row.names=F, col.names=F, sep="\t", na="NA")

q()
n
