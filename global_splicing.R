library(ggplot2)

library(Hmisc)
library(corrplot)
library(plyr)
library(dplyr)

## table with frequencies of alternative splicing and aberrant splicing totals
tab <- read.table("absplicing_asplicing_expr_ggplot.txt",sep="\t",header=TRUE);
tab <- read.table("/Users/naqvia/Desktop/DIPG/test.txt",sep="\t",header=TRUE);

tab <- read.table("/Users/naqvia/Desktop/DIPG/absplicing_asplicing_expr_ggplot.v2.txt",,sep="\t",header=TRUE)

data_long <- gather(tab, Splicing, Percentage, Alternative:Aberrant, factor_key=TRUE)

data_long$Sample <- factor(data_long$Sample, levels =rev(unique(data_long$Sample[order(data_long$Percentage)])),ordered=TRUE)

# Change the colors manually
p <- ggplot(data=data_long, aes(x=Sample, y=Percentage, fill=Splicing)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  scale_fill_manual(values=c('lightblue','blue')) +
  theme_Publication()





##input tables
tab <- read.csv("/Users/naqvia/Desktop/DIPG/tab_global_splicing.txt",sep=",", header=TRUE)
tab_total <-read.csv("/Users/naqvia/Desktop/DIPG/tab_global_splicing.total.txt", sep=",",header=TRUE)

tab_total <- read.csv("/Users/naqvia/Desktop/DIPG/tab_global_splicing.total_updated.txt",sep=",",header=TRUE)
# #historgram of splicing changes in DIPG samples
p<-ggplot(tab_total, aes(x=Counts))  + geom_histogram(color="black", fill="blue",binwidth=1) +  ylab("Differential LSV count ") + xlab("Num of samples containing LSV") + 
  theme_Publication()

# Add mean line
#p+ geom_vline(aes(xintercept=mean(Counts)),
#            color="blue", linetype="dashed", size=1)

p<-ggplot(tab, aes(x=Counts))+ ylab("Differential LSV count ") + xlab("Num of comparisons that agree") +
  geom_histogram(color="black", fill="red",binwidth=2)+
  facet_grid(Subtype ~ .)

mu <- ddply(tab, "Subtype", summarise, grp.mean=mean(Counts))
p+geom_vline(data=mu, aes(xintercept=grp.mean, color="mean"),
             linetype="dashed")

## ven between normal aberrant splicing vs dmg
## venn diagrams
library(VennDiagram)
venn_data <- read.table("/Users/naqvia/Downloads/normal_vs_normal_ped/total_lsvs_norms_dipg.txt", sep="\t",header=TRUE)

venn.diagram(list("Normals" = venn_data$Normals, "DMG" = venn_data$DMG) ,fill = c("lightblue", "darkblue"),
             alpha = c(0.8, 0.8), cex = 2,cat.fontface = .5,lty =2, 
             total.population = 18546,
             main = "Normals vs DMG",
             filename = "/Users/naqvia/Desktop/normal_dmg_venn.tiff")


## correlation between universal mis-spliced genes and their associated expression
tab = read.csv("/Users/naqvia/Desktop/DIPG/test.tmp.txt",header=TRUE)
tab = read.table("/Users/naqvia/Desktop/DIPG/test.tmp.2.txt",header=TRUE, sep="\t")

tab = read.csv("/Users/naqvia/Desktop/DIPG/max_dpsi_tab_for_univ.csv", header=TRUE,sep=",")

## remove columns that have stdv of 0
tab <- tab[, sapply(tab, function(x) { sd(x) != 0} )]

res <- cor(tab)
#corrplot(res, type = "upper", order = "hclust", tl.col = "black", tl.cex = 0.5,tl.srt = 100)
res2<-rcorr(as.matrix(tab))

## Plot corr. matrix with insignificant correlations are leaved blank
corrplot(res2$r, type="lower", order="hclust", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank", tl.col = "black", tl.cex = 0.5,tl.srt = 100)

## table of significant correlations
tab_sign_outlier = read.table("/Users/naqvia/Desktop/DIPG/uni_dpsi_pvals_corrs.sign_outlier.txt",sep="\t",header = TRUE)
head(tab_sign_outlier)

write.table(res2$P, file = "/Users/naqvia/Desktop/DIPG/pval_univ.txt", sep = "\t",row.names = TRUE, col.names = TRUE)
write.table(res2$r, file = "/Users/naqvia/Desktop/DIPG/rcorr_univ.txt", sep = "\t",row.names = TRUE, col.names = TRUE)

## venn diagrams
library(VennDiagram)
venn_data <- read.table("/Users/naqvia/Desktop/DIPG/lsv_list.global.both.txt", sep="\t",header=TRUE)

#venn.diagram(list("Non_H3K28" = venn_data$Non_H3K28, "H3K28" = venn_data$H3K28) ,fill = c("violet", "red"),
#             alpha = c(0.8, 0.8), cex = 2,cat.fontface = .5,lty =2, 
#             total.population = 340,
#             main = "Non_H3K28 vs H3K28",
#             filename = "/Users/naqvia/Desktop/DIPG_venn.tiff")


venn_data <- read.table("/Users/naqvia/Desktop/DIPG/lsv_list.2more.both.txt", sep="\t",header=TRUE)

#venn.diagram(list("Non_H3K28" = venn_data$Non_H3K28, "H3K28" = venn_data$H3K28) ,fill = c("violet", "red"),
#             alpha = c(0.8, 0.8), cex = 2,cat.fontface = .5,lty =2, 
#             main = "Non_H3K28 vs H3K28",
#             filename = "/Users/naqvia/Desktop/DIPG_venn_rec.tiff")

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
