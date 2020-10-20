if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("bladderbatch")

library(sva)

library(bladderbatch)
data(bladderdata)
dat <- bladderEset[1:50,]

pheno = pData(dat)
edata = exprs(dat)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer), data=pheno)

# parametric adjustment
combat_edata1 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

# non-parametric adjustment, mean-only version
combat_edata2 = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=FALSE, mean.only=TRUE)

# reference-batch version, with covariates
combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)

batch <- combined_clin$batch # batches in your case polyA and stranded
corrected_mat <- ComBat_seq(counts = as.matrix(combined_mat), 
                            batch = batch, 
                            group = NULL) # optional, if you have any groups, use them here as a factor

##check to see affect of stranded vs unstranded
gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/tpm_ordered_by_polyA_total.csv",row.names=1, header=TRUE)
batch <- c(rep(1, 52), rep(2, 123))

adjusted <- ComBat_seq(as.matrix(log2(gene_tpms+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted) 
x=plotMDS(log(corrected_mat), cex.lab= 1, cex = 1, col = c( c(rep("blue",52)), c(rep("red",123))), main = "PCA of polyA vs non_polyA", pch=c(rep(18,182)))

x=plotMDS(log(corrected_mat), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("darkred",28 )),c(rep("blue",101)), c(rep("red",6)), c(rep("darkred", 16))), main = "PCA of HGGs vs DMGs", pch=c(rep(16,52), c(rep(17,123))))
legend("bottomright", legend=c("HGGs-H3 WT PolyA","DMG-H3 WT PolyA","DMG-H3 K28 PolyA","HGGs-H3 WT Non-PolyA","DMG-H3 WT Non-PolyA","DMG-H3 K28 Non-PolyA" ),col = c("blue","red","darkred","blue","red","darkred"), pch = c(16,16,16,17,17,17), horiz=TRUE, cex=0.5)



##use 
gene_psi <- read.csv("/Users/naqvia/Desktop/DIPG/psi_ijc_ordered_by_polyA_total.freq36.csv", row.names=1, header=TRUE)
batch_psi<- c(rep(1, 52), rep(2, 123)) ## polyA vs non-polyA
adjusted_psi <- ComBat_seq(as.matrix((gene_psi)), batch=batch_psi, group=NULL)
corrected_mat_psi <- (adjusted_psi) 
x=plotMDS(log(corrected_mat_psi), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("darkred",28 )),c(rep("blue",101)), c(rep("red",6)), c(rep("darkred", 16))), main = "PCA of HGG vs DMG IJC", pch=c(rep(16,52), c(rep(17,123))))
legend("bottomright", legend=c("HGGs-H3 WT PolyA","DMG-H3 WT PolyA","DMG-H3 K28 PolyA","HGGs-H3 WT Non-PolyA","DMG-H3 WT Non-PolyA","DMG-H3 K28 Non-PolyA" ),col = c("blue","red","darkred","blue","red","darkred"), pch = c(16,16,16,17,17,17),horiz=TRUE, cex=0.5)

#gene_psi <- read.csv("/Users/naqvia/Desktop/DIPG/psi_ordered_by_polyA_total.freq36.csv", row.names=1, header=TRUE)
#adjusted_psi <- ComBat_seq(as.matrix((gene_psi)), batch=batch_psi, group=NULL)
#corrected_mat_psi <- (adjusted_psi) 
#adjusted_psi <- ComBat_seq(as.matrix(log2(gene_psi+1)), batch=batch_psi, group=NULL)
#corrected_mat_psi <- 2^(adjusted_psi)
#x=plotMDS(log(corrected_mat_psi), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("darkred",28 )),c(rep("blue",101)), c(rep("red",6)), c(rep("darkred", 16))), main = "PCA of HGGs vs DMGs", pch=c(rep(16,52), c(rep(17,123))), gene.selection = "pairwise")
#legend("bottomright", legend=c("HGGs-H3 WT PolyA","DMG-H3 WT PolyA","DMG-H3 K28 PolyA","HGGs-H3 WT Non-PolyA","DMG-H3 WT Non-PolyA","DMG-H3 K28 Non-PolyA" ),col = c("blue","red","darkred","blue","red","darkred"), pch = c(16,16,16,17,17,17),horiz=TRUE, cex=0.5)


x=plotMDS(log(corrected_mat_psi), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("red",28 )), c(rep("blue",101)), c(rep("red",6)), c(rep("red", 16))), main = "PCA of HGGs (non-DMGs) vs DMGs", pch=c(rep(16,52), c(rep(17,123))), gene.selection = "pairwise")
legend("bottomright", legend=c("HGGs-H3 WT PolyA","DMG-H3 WT PolyA","DMG-H3 K28 PolyA","HGGs-H3 WT Non-PolyA","DMG-H3 WT Non-PolyA","DMG-H3 K28 Non-PolyA" ),col = c("blue","red","red","blue","red","red"), pch = c(16,16,16,17,17,17),horiz=TRUE, cex=0.5)






x=plotMDS((corrected_mat_psi), cex.lab= 1, cex = 1, col = c( c(rep("black", 19)),c(rep("blue",5)), c(rep("purple",28 )), c(rep("red",101)), c(rep("darkred",6)), c(rep("darkgreen", 16))), main = "PCA of polyA vs non-polyA", pch=c(rep(18,175)), top=100)



gene_psi <- read.csv("/Users/naqvia/Desktop/DIPG/psi_ijc_ordered_by_polyA_total.freq2.csv", row.names=1, header=TRUE)
batch_psi<- c(rep(1, 52), rep(2, 123))

adjusted_psi <- ComBat_seq(as.matrix((gene_psi)), batch=batch_psi, group=NULL)
corrected_mat_psi <- (adjusted_psi) 

#adjusted_psi <- ComBat_seq(as.matrix(log2(gene_psi+1)), batch=batch_psi, group=NULL)
#corrected_mat_psi <- 2^(adjusted_psi) 


x=plotMDS(log(corrected_mat_psi), cex.lab= 1, cex = 1, col = c( c(rep("blue",52)), c(rep("red",123))), main = "PCA of polyA vs non-polyA", pch=c(rep(18,175)))
x=plotMDS(log(corrected_mat_psi), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("red",28 )), c(rep("blue",101)), c(rep("red",6)), c(rep("red", 16))), main = "PCA of HGGs vs DMGs", pch=c(rep(18,52), c(rep(19,123))))

polyA_only <- (gene_psi[,1:53])
nonpolyA   <- (gene_psi[,54:123])
x=plotMDS(log(nonpolyA), cex.lab= 1, cex = 1, col = c( c(rep("blue",19)), c(rep("red",33))), main = "PCA of HGGs vs DMGs (polyA)", pch=c(rep(18,52)))



polyA_only

