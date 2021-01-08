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

## DIPG table with unipro results
file="/Users/naqvia/Desktop/DIPG/output/dominant_events_lsvs.total.intersectUnip.ggplot.txt" 
file="/Users/naqvia/Desktop/DIPG/re_analysis/dominant_events_lsvs.total.intersectUnip.uniq.ggplot.txt" 

file="/Users/naqvia/Desktop/DIPG/re_analysis/tables/dominant_events_lsvs.total.intersectUnip_and_neg.uniq.ggplot.v3.txt"

dpsi_unip <- read.table(file, header=TRUE,sep = "\t")

set.seed(123)

ggstatsplot::ggbetweenstats(
  data = dpsi_unip, 
  x = Uniprot, 
  y = dPSI,
  k = 2,
  notch = FALSE,
  mean.ci = TRUE,
  outlier.tagging = TRUE,
  outlier.label = LSV,
  title = "Distribution of dPSI across functional sites",
  #type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = TRUE + 
  geom_text_repel()
) + theme_Publication()

# 
file="/Users/naqvia/Desktop/DIPG/re_analysis/dominant_events_lsvs.total.intersectUnip.ts.uniq.ggplot.txt" 
dpsi_unip_ts <- read.table(file, header=TRUE,sep = "\t")

set.seed(123)

plot_ts <- ggstatsplot::ggbetweenstats(
  data = dpsi_unip_ts, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  #outlier.label = LSV, # label to attach to outlier values
  #outlier.label.args = list(color = "red"), # outlier point label color
  notch = FALSE,
  mean.ci = TRUE,
  outlier.tagging = FALSE,
  title = "Tumor supressors",
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) 

file="/Users/naqvia/Desktop/DIPG/re_analysis/dominant_events_lsvs.total.intersectUnip.onco.uniq.ggplot.txt" 
dpsi_unip_onco <- read.table(file, header=TRUE,sep = "\t")

set.seed(123)

plot_onco <- ggstatsplot::ggbetweenstats(
  data = dpsi_unip_onco, 
  x = Uniprot, 
  y = dPSI,
  k = 3,
  nboot = 15,
  #outlier.label = LSV, # label to attach to outlier values
  #outlier.label.args = list(color = "red"), # outlier point label color
  notch = FALSE,
  mean.ci = TRUE,
  outlier.tagging = FALSE,
  title = "Oncogenes",
  type = "robust",
  xlab = "Unipro-defined Site",
  pairwise.comparisons = TRUE,
  messages = FALSE
) 


ggstatsplot::combine_plots(
  plot_ts, plot_onco,
  nrow = 2,
  title.text = "Distribution of dPSI across functional sites in cancer genes ",
  title.size = 14,
  caption.size = 12
)

## CLK AS events (rMATS)
file="/Users/naqvia/Desktop/DIPG/CLK_analysis/top_analysis/tab_for_ggplot.neg.v2.txt"

dpsi_unip <- read.table(file, header=TRUE,sep = "\t")

set.seed(123)

ggstatsplot::ggbetweenstats(
  data = dpsi_unip, 
  x = Type, 
  y = dPSI,
  k = 2,
  notch = FALSE,
  mean.ci = FALSE,
  outlier.tagging = TRUE,
  outlier.label = LSV,
  title = "Distribution of dPSI across functional sites",
  #type = "robust",
  xlab = "Type",
  pairwise.comparisons = FALSE,
  messages = TRUE 
)


file="/Users/naqvia/Desktop/DIPG/isoform_switch_pie.dat"
switch_types <- read.table(file, header=TRUE,sep = "\t")

# Basic piechart
ggplot(switch_types, aes(x="", y=Counts, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_minimal() +
  scale_fill_manual(values=c("darkblue", "red")
  #theme_Publication() # remove background, grid, numeric labels
  #scale_fill_brewer(palette="Purples")


# Load ggplot2
library(ggplot2)

# Create Data
data <- data.frame(
  group=LETTERS[1:5],
  value=c(13,7,9,21,2)
)

file="/Users/naqvia/Desktop/DIPG/re_analysis/dominant_events_lsvs.total.intersectUnipMod.piechart_nums.dat" 
mods_types <- read.table(file, header=TRUE,sep = "\t")

# Basic piechart
ggplot(mods_types, aes(x="", y=Counts, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_minimal() +
  # theme_Publication() # remove background, grid, numeric labels
 scale_fill_brewer(palette="Purples")



file="/Users/naqvia/Desktop/DIPG/re_analysis/dominant_events_lsvs.total.intersectUnipMod.onco.piechart_nums.dat" 
mods_types <- read.table(file, header=TRUE,sep = "\t")

# Basic piechart
ggplot(mods_types, aes(x="", y=Counts, fill=Type)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) + theme_minimal() +
  # theme_Publication() # remove background, grid, numeric labels
  scale_fill_brewer(palette="Purples")

file="/Users/naqvia/Desktop/DIPG/isoform_switching_vs_same.dat" 
isoform_switch_nums <- read.table(file, header=TRUE,sep = "\t")

file="/Users/naqvia/Desktop/DIPG/re_analysis/isoform_switching_vs_same_updated.dat"
isoform_switch_nums <- read.table(file, header=TRUE,sep = "\t")

isoform_switch_nums$Sample <- factor(isoform_switch_nums$Sample, levels =rev(unique(isoform_switch_nums$Sample[order(isoform_switch_nums$Counts)])),ordered=TRUE)
#isoform_switch_nums$Sample  # notice the changed order of factor levels

# plot
ggplot(data=isoform_switch_nums, aes(x=Sample, y=Counts, fill=Type)) + 
  geom_bar(stat="identity")  +
  scale_fill_manual(values=c('darkblue','red')) +
  theme_Publication()

#isoform_switching_vs_non.ts.txt
file="/Users/naqvia/Desktop/DIPG/isoform_switching_vs_strong-non.ts.txt" 
isoform_switch_nums <- read.table(file, header=TRUE,sep = "\t")

isoform_switch_nums$Sample <- factor(isoform_switch_nums$Sample, levels =rev(unique(isoform_switch_nums$Sample[order(isoform_switch_nums$Counts)])),ordered=TRUE)
#isoform_switch_nums$Sample  # notice the changed order of factor levels

# plot
ggplot(data=isoform_switch_nums, aes(x=Sample, y=Counts, fill=Type)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('darkblue','red')) +
  theme_Publication()

file="/Users/naqvia/Desktop/DIPG/isoform_switching_vs_non.onco.txt" 
isoform_switch_nums <- read.table(file, header=TRUE,sep = "\t")

isoform_switch_nums$Sample <- factor(isoform_switch_nums$Sample, levels =rev(unique(isoform_switch_nums$Sample[order(isoform_switch_nums$Counts)])),ordered=TRUE)
#isoform_switch_nums$Sample  # notice the changed order of factor levels

# plot
ggplot(data=isoform_switch_nums, aes(x=Sample, y=Counts, fill=Type)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c('darkblue','red')) +
  theme_Publication()

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