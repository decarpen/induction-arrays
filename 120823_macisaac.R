library(gridExtra)
library(ggplot2)
library(reshape2)
library(impute)
library(GOstats)
library(GO.db)
library(org.Sc.sgd.db)
library(GSEABase)
library(xtable)
source("analysis_fxns.R")
source("heatmap_fxn.R")

#p.dat = read.delim("SVD_norm_carbon.txt")
p.dat = read.delim("SVD_norm_phosphate.txt")


mac.dat = scan(file="~/Documents/SCIENCE/BOTSTEIN/EXPERIMENTS/120601_ZEV_SNF1_screen/DATA/MacIsaac_2006_data/orfs_by_factor/orfs_by_factor_none_cons0.txt",
               what=character(0),sep="\n")
mac.dat = strsplit(mac.dat,"\t",fixed=T)

tf.dat = list()
for (i in 1:length(mac.dat)) {
  tf.dat[[ mac.dat[[i]][1] ]] = mac.dat[[i]][2:length(mac.dat[[i]])]
}
adr1 = tf.dat$ADR1
sip4 = tf.dat$SIP4
mig1 = tf.dat$MIG1
rgt1 = tf.dat$RGT1
msn2 = tf.dat$MSN2

print("YPL248C" %in% mig1)

mdat = matrix(F,nrow=nrow(p.dat),ncol=4)
rownames(mdat) = rownames(p.dat)
colnames(mdat) = c("ADR1", "MIG1", "RGT1", "SIP4")
stuff = list(adr1, mig1, rgt1, sip4)
for (i in 1:4) {
  sel = rownames(p.dat) %in% stuff[[i]]
  mdat[,i] = sel
}
mdat = mdat+0
rownames(mdat) = orf2name(rownames(mdat))
g = ggheat("promoter_presence_phosphate.pdf",mdat,c("white","white","cornflower blue"))

