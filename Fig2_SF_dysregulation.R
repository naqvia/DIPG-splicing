library("sva")
library("EnhancedVolcano")
library("DESeq2")
library("ggplot2")
library("Hmisc")
library("corrplot")
library("plyr")
library("dplyr")
library("VennDiagram")

gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/tpm_ordered_by_polyA_dmg_vs_5norms.SFs_only.csv",row.names=1, header=TRUE)
batch <- c(rep(1, 33), rep(2, 22), rep(1,5))

filtered.counts <- countTable[rowSums(corrected_mat>=10) > 10, ]
countTable <- filtered.counts ## filter out  genes with low read counts

## batch correction
adjusted <- ComBat_seq(as.matrix(log2(gene_tpms+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted) 

## convert to dataframe and integer for DESeq2
corrected_df <- as.data.frame(corrected_mat)
corrected_df[] <- sapply( corrected_df, as.integer )

##check to see if batch correction worked?
x=plotMDS(log(corrected_mat), cex.lab= 1, cex = 1, col = c( c(rep("red", 33)),c(rep("darkred",22)),  c(rep("blue",5))), main = "PCA of  DMGs vs Norms Expr", pch=c( rep(16,33), rep(17,22), rep(16,5)))
          
## construct metadata
design = data.frame(row.names = colnames(corrected_df),
                    condition = c(rep("DMG",55), rep("Norm",5)),  
                    libType   = c(rep("paired-end",60)))

singleSamples = design$libType == "paired-end"
new_countTable = (corrected_df[ , singleSamples ])
condition = design$condition[ singleSamples ]

cds = newCountDataSet( (new_countTable), condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )

cds = estimateDispersions(cds) 

res = nbinomTest( cds, "DMG", "Norm")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")

## plot 
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                ylim = c(0,20),
                xlim =c(-5,5),
                title = 'DMG versus Norm',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 2.0)



