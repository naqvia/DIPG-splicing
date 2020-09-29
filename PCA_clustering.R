library(rgl)
library(gplots)
library(reshape2)
library(ggplot2)
library(limma)



## DMGs only
file = "/Users/naqvia/Desktop/DIPG/psi_tumor.dpsi.clk1.csv"
psi_dpsi <- read.csv(file ,header=TRUE,row.names=1)
x=plotMDS((psi_dpsi), cex.lab= 1, cex = 1, col = c( c(rep("blue",9)), c(rep("red",39)) ), main = "PCA of non-CLK vs CLK PSI", pch=c(rep(18,47)))

## HGGs + DMGs only
file = "/Users/naqvia/Desktop/DIPG/psi_tumor.dpsi.clk_wHGGs.csv"
psi_dpsi <- read.csv(file ,header=TRUE,row.names=1)
x=plotMDS((psi_dpsi), cex.lab= 1, cex = 1, col = c( c(rep("blue",39)), c(rep("purple",98)),c(rep("darkgreen",9)), c(rep("red",39)) ), main = "PCA of non-CLK vs CLK PSI in HGGs", pch=c(rep(18,184)))

## PCA of HGGs (non-DMG) without H3K28M (blue) and DMGs WT (red) and DMGs H3K28 (red) based on recurrent splicing
file = "/Users/naqvia/Desktop/DIPG/deltapsi_tab_rMATS.SE.HGGs_and_DMGs.sign05_rec36.recalc_h3k28.csv"

psi_lsv <- read.csv(file ,header=TRUE,row.names=1)
x=plotMDS((psi_lsv), cex.lab= 1, cex = 1, col = c( c(rep("blue",128)), c(rep("red",11)),c(rep("darkred",35)) ), main = "PCA of HGGs PSI", pch=c(rep(18,139),rep(20,35)))
legend("bottomright", legend=c("HGGs-H3 WT","DMG-H3 WT","DMG-H3 K28"),col = c("blue","red","darkred"), pch = c(18,18,20),horiz=TRUE, cex=0.8)

##PCA on expression