library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(ggplot2)
library(dplyr)
library(ggstatsplot)


library(ggplot2)
library(dplyr)
library(ggstatsplot)

## micro vs non-micro
file ="/Users/naqvia/Desktop/DIPG/output/recurrent/dominant_events_lsvs.isoform_switches.non_micro_micro.recurrent10.intersectUnip.uniq.ggplot.txt"
dpsi_unip <- read.table(file, header=TRUE,sep = "\t")

file = "/Users/naqvia/Desktop/DIPG/re_analysis/dominant_events_lsvs_simple_recurrent2.len.micro_vs_cano.ggplot.txt"

dpsi_unip <- read.table(file, header=TRUE,sep = "\t")

library(ggstatsplot)
library(gapminder)

# for reproducibility
set.seed(123)

ggstatsplot::ggbetweenstats(
  data = dpsi_unip, 
  x = Type, 
  y = dPSI,
  k = 2,
  caption = "Mann–Whitney U test",
  xlab = "Exon Type",
  notch = FALSE,
  mean.ci = TRUE,
  outlier.tagging = FALSE,
  title = "Distribution of dPSI across functional sites",
  type="np",
  conf.level = 0.95,
  pairwise.comparisons = TRUE,
  messages = TRUE
) + theme_Publication()

file = "/Users/naqvia/Desktop/DIPG/re_analysis/dominant_events_lsvs.micro.intersectUnip.uniq.ggplot.txt"
psi_unip <- read.table(file, header=TRUE,sep = "\t")
# for reproducibility
set.seed(123)

ggstatsplot::ggbetweenstats(
  data = dpsi_unip, 
  x = Uniprot, 
  y = dPSI,
  k = 2,
  caption = "Mann–Whitney U test",
  xlab = "Uniprot Site",
  notch = FALSE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  outlier.label = LSV,
  title = "Distribution of dPSI across functional sites",
  type="np",
  conf.level = 0.95,
  pairwise.comparisons = TRUE,
  messages = TRUE
) + theme_Publication()

tab <- read.table("/Users/naqvia/Desktop/DIPG/len_me.dat",  header=TRUE,sep="\t")
tab_onco <- read.table("/Users/naqvia/Desktop/DIPG/len_me.onco.dat", header=TRUE,sep="\t")
tab_ts <- read.table("/Users/naqvia/Desktop/DIPG/len_me.ts.dat", header=TRUE,sep="\t")

ggplot(tab, aes(Length,LSV.Counts)) + 
  geom_bar(stat="identity", fill="blue") + xlab("Exon length") + ylab("LSV Counts") + 
  geom_text(aes(label=LSV.Counts), vjust=-0.2, size=2.5)  +
  theme_Publication()


geom_text(aes(label=len), vjust=1.6, color="white", size=3.5)+
  theme_minimal()


,type = "b", frame = FALSE, pch = 18, col = "darkgreen", xlab = "Microexon Length", ylab = "Differential LSV Count") 
plot(tab_onco$Length,tab_onco$LSV.Counts,type = "b", frame = FALSE, pch = 18, col = "red", lty = 1, xlab = "Microexon Length", ylab = "Differential LSV Count") 

# Add a second line
lines(tab_onco$Length,tab_onco$LSV.Counts, pch = 18, col = "red", type = "b", lty = 2)

lines(tab_ts$Length,tab_ts$LSV.Counts, pch = 19, col = "blue", type = "b", lty = 1)

legend("topleft", legend=c("Oncogenes", "Tumor suppressors"),
       col=c("red", "blue"), lty = 1:1, cex=0.8)


frame_total.dat
frame_total.onco.dat
frame_total.ts.dat

file="/Users/naqvia/Desktop/DIPG/frame_total.dat"
frame_total <- read.table(file, header=TRUE,sep = "\t")

file_onco="/Users/naqvia/Desktop/DIPG/frame_total.onco.dat"
frame_onco <- read.table(file_onco, header=TRUE,sep = "\t")

file_ts="/Users/naqvia/Desktop/DIPG/frame_total.ts.dat"
frame_ts <- read.table(file_ts, header=TRUE,sep = "\t")

# Basic piechart
ggplot(frame_total, aes(x="", y=Counts, fill=Frame)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_minimal() +
  scale_fill_manual(values=c("blue", "red"))

ggplot(frame_onco, aes(x="", y=Counts, fill=Frame)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_minimal() +
  scale_fill_manual(values=c("blue", "red"))

ggplot(frame_ts, aes(x="", y=Counts, fill=Frame)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_minimal() +
  scale_fill_manual(values=c("blue", "red"))

#theme_Publication() # remove background, grid, numeric labels
#scale_fill_brewer(palette="Purples")

file="/Users/naqvia/Desktop/DIPG/frame_me_len_ggplot.dat"
me_frame <- read.csv(file, header=TRUE,sep = ",")

me_frame$type <- factor(me_frame$type, levels =rev(unique(me_frame$type[order(me_frame$counts)])),ordered=TRUE)
#isoform_switch_nums$Sample  # notice the changed order of factor levels

# plot
ggplot(data=me_frame, aes(x=type, y=counts, fill=frame)) + 
  geom_bar(stat="identity")  +
  scale_fill_manual(values=c('blue','red')) +
  theme_Publication()


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

