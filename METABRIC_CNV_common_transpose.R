# To get transpose of METABRIC CNV table
cat("\014")  
rm(list = ls())
setwd("F:/BioInformatics/Paper Data")
mydata<- read.table("METABRIC_CNV_common.txt")
class(mydata)
str(mydata)
summary(mydata)
head(mydata)
mydata.t = t(mydata)
ncol(mydata)
nrow(mydata)
rownames(mydata)
colnames(mydata)
mydata[1,]
mydata[,1]
write.table(mydata.t,"METABRIC_CNV_common_transpose.txt",quote = FALSE, row.names = FALSE, col.names = FALSE)
 