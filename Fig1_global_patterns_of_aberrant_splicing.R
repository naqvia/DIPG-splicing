library("sva")
library("EnhancedVolcano")
library("DESeq2")
library("ggplot2")
library("Hmisc")
library("corrplot")
library("plyr")
library("dplyr")
library("VennDiagram")

## histogram of differential LSVs in cohort
tab_total <- read.csv("/Users/naqvia/Desktop/DIPG/tab_global_splicing.total_updated.txt",sep=",",header=TRUE)
p<-ggplot(tab_total, aes(x=Counts))  + geom_histogram(color="black", fill="blue",binwidth=1) +  ylab("Differential LSV count ") + xlab("Num of samples containing LSV") + 
  theme_Publication()

## overlap normals vs DMG LSVs
venn_data <- read.table("/Users/naqvia/Downloads/normal_vs_normal_ped/total_lsvs_norms_dipg.txt", sep="\t",header=TRUE)
venn.diagram(list("Normals" = venn_data$Normals, "DMG" = venn_data$DMG) ,fill = c("lightblue", "darkblue"),
             alpha = c(0.8, 0.8), cex = 2,cat.fontface = .5,lty =2, 
             total.population = 18546,
             main = "Normals vs DMG",
             filename = "/Users/naqvia/Desktop/normal_dmg_venn.tiff")


## generate correlation matrix of universally mis-spliced events
tab = read.csv("/Users/naqvia/Desktop/DIPG/max_dpsi_tab_for_univ.csv", header=TRUE,sep=",")

## remove columns that have stdv of 0
tab <- tab[, sapply(tab, function(x) { sd(x) != 0} )]

res2<-rcorr(as.matrix(tab))

## Plot corr. matrix with insignificant correlations are leaved blank
p <- corrplot(res2$r, type="lower", order="hclust", 
         p.mat = res2$P, sig.level = 0.05, insig = "blank", tl.col = "black", tl.cex = 0.7,tl.srt = 1) 

# write table txt of significant correlations
# write.table(res2$P, file = "/Users/naqvia/Desktop/DIPG/pval_univ.txt", sep = "\t",row.names = TRUE, col.names = TRUE)
# write.table(res2$r, file = "/Users/naqvia/Desktop/DIPG/rcorr_univ.txt", sep = "\t",row.names = TRUE, col.names = TRUE)


## cluster (PCA) with batch correction
##based on gene expression
gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/tpm_ordered_by_polyA_total.csv",row.names=1, header=TRUE)

## batch correction
batch <- c(rep(1, 52), rep(2, 123))
adjusted <- ComBat_seq(as.matrix(log2(gene_tpms+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted) 

x=plotMDS(log(corrected_mat), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("darkred",28 )),c(rep("blue",101)), c(rep("red",6)), c(rep("darkred", 16))), main = "PCA of HGGs vs DMGs Expr", pch=c(rep(16,52), c(rep(17,123))))
legend("bottomright", legend=c("HGGs-H3 WT PolyA","DMG-H3 WT PolyA","DMG-H3 K28 PolyA","HGGs-H3 WT Non-PolyA","DMG-H3 WT Non-PolyA","DMG-H3 K28 Non-PolyA" ),col = c("blue","red","darkred","blue","red","darkred"), pch = c(16,16,16,17,17,17), horiz=TRUE, cex=0.5)

##based on splicing
gene_psi <- read.csv("/Users/naqvia/Desktop/DIPG/psi_ijc_ordered_by_polyA_total.freq36.csv", row.names=1, header=TRUE)
batch_psi<- c(rep(1, 52), rep(2, 123)) ## polyA vs non-polyA

adjusted_psi <- ComBat_seq(as.matrix((gene_psi)), batch=batch_psi, group=NULL)
corrected_mat_psi <- (adjusted_psi) 

x=plotMDS(log(corrected_mat_psi), cex.lab= 1, cex = 1, col = c( c(rep("blue", 19)),c(rep("red",5)), c(rep("darkred",28 )),c(rep("blue",101)), c(rep("red",6)), c(rep("darkred", 16))), main = "PCA of HGG vs DMG IJC", pch=c(rep(16,52), c(rep(17,123))))
legend("bottomright", legend=c("HGGs-H3 WT PolyA","DMG-H3 WT PolyA","DMG-H3 K28 PolyA","HGGs-H3 WT Non-PolyA","DMG-H3 WT Non-PolyA","DMG-H3 K28 Non-PolyA" ),col = c("blue","red","darkred","blue","red","darkred"), pch = c(16,16,16,17,17,17),horiz=TRUE, cex=0.5)

## make volcano plot of WT-H3 vs H3K28M for differential expression analysis
gene_tpms <- read.csv("/Users/naqvia/Desktop/DIPG/tpm_ordered_by_polyA_dmg_only.csv",row.names=1, header=TRUE)
batch <- c(rep(1, 33), rep(2, 22))
filtered.counts <- countTable[rowSums(corrected_mat>=10) > 10, ]
countTable <- filtered.counts ## filter out  genes with low read counts

## batch correction
adjusted <- ComBat_seq(as.matrix(log2(gene_tpms+1)), batch=batch, group=NULL)
corrected_mat <- 2^(adjusted) 

## convert to dataframe and integer for DESeq2
corrected_df <- as.data.frame(corrected_mat)
corrected_df[] <- sapply( corrected_df, as.integer )

## construct metadata
design = data.frame(row.names = colnames(corrected_df),
                    condition = c(rep("H3-WT",5), rep("H3K28",28), rep("H3-WT", 6), rep("H3K28", 16) ),  
                    libType   = c(rep("paired-end",55)))

singleSamples = design$libType == "paired-end"
new_countTable = (corrected_df[ , singleSamples ])
condition = design$condition[ singleSamples ]

cds = newCountDataSet( (new_countTable), condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )

cds = estimateDispersions(cds) 

res = nbinomTest( cds, "H3-WT", "H3K28")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")

## plot 
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                ylim = c(0,2.5),
                xlim =c(-3,3),
                title = 'H3-WT versus H3K27',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 2.0)





##theme for all plots
theme_Publication <- function(base_size=15, base_family="Helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            legend.direction = "vertical",
            legend.key.size= unit(0.5, "cm"),
            # legend.margin = unit(0.5, "cm"),
            legend.margin = margin(5,5,5,5),
            legend.title = element_text(face="bold"),
            #plot.margin=unit(c(10,5,5,5),"mm"),
            plot.margin=unit(c(10,5,5,10),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}