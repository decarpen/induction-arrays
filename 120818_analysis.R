library(ggplot2)
library(reshape2)
library(impute)
source("analysis_fxns.R")
source("heatmap_fxn.R")

rawdata = read.delim("~/Documents/SCIENCE/BOTSTEIN/EXPERIMENTS/120601_ZEV_SNF1_screen/DATA/120817_PUMA_dld.pcl")

#drop gene weights
rawdata = rawdata[-1,-3]

#remove row names
data = rawdata[!duplicated(rawdata[,2]),]
rownames(data) = data[,2]
data = as.matrix(data[,3:ncol(data)])

#adjust column names
strains = rep(c("12416", "12424", "12416", "12424"), each=8)
media = rep(c("P","C"), each=16)
time = rep(c("0","5","10","15","20","30","45","60"),4)
clist = paste(strains,media,time,sep="_")
colnames(data) = clist
rm(strains,media,time,clist)

#impute missing data
#imp.dat = impute.knn(data)$data
#write.table(imp.dat, "imputed_data.txt", quote=F, sep="\t")
imp.dat = read.delim("imputed_data.txt")

#zero-normalize data
p.12416 = imp.dat[,1:8]
p.12424 = imp.dat[,9:16]
c.12416 = imp.dat[,17:24]
c.12424 = imp.dat[,25:32]
zero.dat = cbind(zeroNorm(p.12424), zeroNorm(p.12416), zeroNorm(c.12424), zeroNorm(c.12416))
write.table(zero.dat, "zero_normalized.txt", sep="\t", quote=F, row.names=F)

#remove those genes with less than 2-fold change
diff.dat = dropLow(zero.dat,1)
write.table(diff.dat, "twofoldchange_data.txt", quote=F, sep="\t")

#just each limitation change
diff.p.dat = dropLow(cbind(zeroNorm(p.12424),zeroNorm(p.12416)),1)
diff.c.dat = dropLow(cbind(zeroNorm(c.12424),zeroNorm(c.12416)),1)

p = ggheat("p_lim_twofold.pdf", diff.p.dat, c("blue","white","red"))
p = ggheat("c_lim_twofold.pdf", diff.c.dat, c("blue","white","red"))

#looking at just the zero timepoints.
zeros = imp.dat[,c(9,1,25,17)]
test = apply(zeros,1,function(x) {return(sd(x)>0.5779777)})
p = ggheat("zero_timepoints.pdf", zeros[test,], c("blue","white","red"),width.adjust=2.5)
