library("car")

tab_splicing<-read.table("/Users/naqvia/Desktop/DIPG/dominant_changes_lsv.isoform_switch.recurrent10.druggable.tab_format.uniq.txt",header=TRUE, sep="\t")
scatterplot(tab_splicing$Freq  ~ tab_splicing$dPSI.avg.,data=tab_splicing,
            smoother = FALSE,  frame = FALSE,
            ,col="red")

# library
library(ggplot2)

# Keep 30 first rows in the mtcars natively available dataset
tab_splicing<-read.table("/Users/naqvia/Desktop/DIPG/dominant_changes_lsv.all.druggable.tab_format.uniq.txt",header=TRUE, sep="\t")
data = tab_splicing

# 1/ add text with geom_text, use nudge to nudge the text
ggplot(data, aes(x=dPSI.avg., y=Freq)) +
  geom_point(col="red") + # Show dots
  geom_label(
    label=data$Gene, 
    nudge_x = 0.015, nudge_y = 0.015,cex=2, 
  ) + theme_Publication()

ggplot(data, aes(x=dPSI.avg., y=Freq)) +
  geom_point(col="red") + # Show dots
  geom_label(
    label=data$Gene, 
    nudge_x = 0.015, nudge_y = 0.015,cex=2, 
  ) + theme_Publication()

  geom_text(
    label=data$Gene, 
    nudge_x = 0.25, nudge_y = 0.25, 
    check_overlap = T
  )

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
  