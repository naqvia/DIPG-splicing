library("reshape")
splice_subtype <- read.csv("/Users/naqvia/Desktop/DIPG/splice_freqs_by_subtype.IJC.w_perc.v2.csv", row.names = 1, header=TRUE)




ggplot(splice_subtype, aes(x=H3K28, y=WT)) + geom_count(shape=18) + geom_point(shape=18,color="darkgreen") + 
  ggtitle("Aberrant Splicing Frequencies (%) across WT vs H3K28")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_Publication()

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

BiocManager::install('EnhancedVolcano')

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("airway")

library(airway)
library(magrittr)





full_path = "/Users/naqvia/Desktop/DIPG/gene_counts_tpm.by_status.v2.txt"
countTable <- read.csv(full_path,header=TRUE,row.names=1) 
head(countTable)

filtered.counts <- countTable[rowSums(countTable>=100) > 10, ]
countTable <- filtered.counts 

##construct metadata
design = data.frame(row.names = colnames( countTable ),
                    condition = c(rep("WT",11), rep("H3K27",35) ),  
                    libType   = c(rep("paired-end",46)))

singleSamples = design$libType == "paired-end"
new_countTable = countTable[ , singleSamples ]
condition = design$condition[ singleSamples ]

cds = newCountDataSet( new_countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cds = estimateDispersions( cds,  fitType = c("local"))

res = nbinomTest( cds, "WT", "H3K27")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")


res <- results(cds,
               contrast = c("WT", "H3K27"))
res <- lfcShrink(cds,
                 contrast = c("WT", "H3K27"), res=res, type = 'normal')

EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'pval')
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                ylim = c(0,2.5),
                title = 'WT versus H3K27',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 2.0)

##aes(x = log2FoldChange, y = -log10(pval)))
library("EnhancedVolcano")
library("DESeq2")
full_path = "/Users/naqvia/Desktop/DIPG/gene_counts_tpm.sfactors.rounded_ceil.v2.tsv"
countTable <- read.table(full_path,header=TRUE,row.names=1) 
head(countTable)

##construct metadata
design = data.frame(row.names = colnames( countTable ),
                    condition = c(rep("norm",5), rep("DIPG",48) ),  
                    libType   = c(rep("paired-end",53)))

singleSamples = design$libType == "paired-end"
new_countTable = countTable[ , singleSamples ]
condition = design$condition[ singleSamples ]

cds = newCountDataSet( new_countTable, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
cds = estimateDispersions( cds )

res = nbinomTest( cds, "norm", "DIPG")
res$Significant <- ifelse(res$pval< 0.05, "P-val < 0.05", "Not Sig")

ggplot(res, aes(x = log2FoldChange, y = -log10(pval))) +
  geom_point(aes(color = Significant)) +
  scale_color_manual(values = c("grey", "red")) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom") +
  geom_text_repel(
    data = subset(res, pval <= 0.05 & abs(res$log2FoldChange)>=1),
    aes(label = id),
    size = 2,
    box.padding = unit(.5, "lines"),
    point.padding = unit(.5, "lines"))

##enhanced
EnhancedVolcano(res,
                lab = (res$id),
                x = 'log2FoldChange',
                y = 'pval',
                xlim = c(-5,5),
                ylim = c(0,14),
                title = 'BS vs DMG (SFs)',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 1.0,
                labSize = 3.0)


##theme
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