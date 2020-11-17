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
  
  library(UpSetR)
  library(ggplot2)
  library(gplots)
  
  
  tab_total <- read.csv("/Users/naqvia/Desktop/DIPG/SURVIV_Result_P.sign.anno.dmg_only.v3.txt",sep=",",header=TRUE)
  
  p<-ggplot(tab_total, aes(x=Type))  + geom_bar(color="black", fill="red", stat = "count") +  ylab("SURVIV Exon Count") + xlab("Class") +
    theme_Publication() +   scale_x_discrete(name = " ", limits = c("Other","ONC","TF","EPI","TS",
                                                                    "Microexon","SF","Kinase","SWISWF"))
  
  p<-ggplot(tab_total, aes(x=SpliceID, fill=Type))  + geom_histogram(color="black", fill="blue",binwidth=1) +  ylab("LSV count ") + xlab("Num of samples with Differential LSV")+
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
  
  
  ## upsetR plots for recurrent LSVs
  #lsvs_ball = read.table("/Users/naqvia/Desktop/TARGET/target_dominant_events.recurrent177.coords.uniq.txt",sep="\t")
  #lsvs_lgg = read.table("/Users/naqvia/Desktop/LGG/lgg_dominant_changes.recurrent40.coords.uniq.txt",sep="\t")
  #lsvs_dipg =read.table("/Users/naqvia/Desktop/DIPG/output/dominant_changes_lsv_list.recurrent10.coords.uniq.txt",sep="\t")
  #lsvs_hgg  =read.table("/Users/naqvia/Desktop/HGG/psi_majiq/hgg_dominant_events.recurrent27.coords.uniq.txt", sep="\t")
  
  ## lsvs overlapping protein domains + recurrent
  ball = read.table("/Users/naqvia/Desktop/TARGET/dominant_events_lsvs.ball.recurrent177.intersectUnipDomain.coords.uniq.txt",sep="\t")
  lgg  = read.table("/Users/naqvia/Desktop/LGG/dominant_events_lsvs.lgg.recurrent40.intersectUnipDomain.coords.uniq.txt",sep="\t")
  hgg  =read.table("/Users/naqvia/Desktop/HGG/dominant_events_lsvs.hgg.recurrent27.intersectUnipDomain.coords.uniq.txt", sep="\t")
  dipg =read.table("/Users/naqvia/Desktop/DIPG/dominant_events_lsvs.dipg.recurrent10.intersectUnipDomain.coords.uniq.txt",sep="\t")
  
  
  listInput <- list("B-ALL" =ball$V1, 
                    "LGG"   =lgg$V1,
                    "DIPG"  =dipg$V1,
                    "HGG"   =hgg$V1)
  
  k = upset(fromList(listInput), 
            mainbar.y.label = "LSV Intersections", sets.x.label = "# of Differential LSVs",
            mb.ratio = c(0.5,0.50),text.scale = c(1.3, 1.3, 1.3, 1.3, 2, 1.4),point.size = 3.5, line.size = 2)
  