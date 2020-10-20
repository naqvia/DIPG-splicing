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

