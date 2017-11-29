setwd("C:/Users/luigi/Dropbox/LoReAn_paper/Data/files-paper-LoReAn/V_dahliae/intron/")
library(reshape2)
library(ggplot2)
library(gridExtra)
library(agricolae)
library(ggrepel)

##############################################
######## Introns #############################
JR2_intron <- read.csv("./onlyIntrons.csv", header=T, fill=T)


acclines2 <- data.frame(
  x = c(120, 160, 180),
  y = c(0, 0, 0),
  xend = c(0, 0, 0),
  yend = c(120, 160, 180),
  color = c("black", "red", "orange")
) 

intron.jr <- ggplot(JR2_intron, aes(x=spe_intro, y=sens_intro)) +
  geom_point(aes(shape = masking, color = masking), size =1.25) + 
  #geom_segment(aes(x=x, y=y, xend=xend, yend=yend, color=color), data=acclines2) +
  coord_cartesian(ylim=c(60, 85), xlim=c(60, 80)) +
  geom_label_repel(
    aes(fill = JR2_intron$pipeline, label= ifelse(acc_intro > 75 | acc_intro < 68 | annot == "GeneMark-ES-F" | annot == "VDAG_Jr2_Annotation" | annot == "CodingQuarry", as.character(JR2_intron$annot),"")),
    fontface = 'plain',  color ="black", size= 2.0,
    box.padding= unit(0.35, "lines"),
    point.padding= unit(0.5, "lines"),
    force = 2, 
    segment.color= 'grey50') +
  theme_bw() +
  theme(legend.position = 'none') 

ggsave(file="jr2.introns.pdf", intron.jr, width = 120, height = 120, units = "mm")



###################################
### Intron Analysis  ##############



int_gene <- read.csv2("./Summary_gene_intron.csv", header=T,
                      sep= ",", dec = ".", stringsAsFactors = F)
#int_gene$value <- paste (int_gene$masking, "_", int_gene$pipeline)
int_gene$accuracy1 <- int_gene$accuracy * 100



with(int_gene, cor.test(as.numeric(acc_intro), as.numeric(accuracy1), method=c("spearman")))


int_gene2 <- ggplot(int_gene, aes(x=acc_intro, y=accuracy1)) +
  geom_point(aes(shape = as.character(int_gene$shape), color = masking), size =1.25) +
  #geom_segment(aes(x=x, y=y, xend=xend, yend=yend, color=color), linetype="dashed", data=acclines_r) +
  coord_cartesian(ylim=c(25, 70), xlim=c(60, 80)) +
  geom_label_repel(
    aes(fill = int_gene$pipeline, label= ifelse(acc_intro > 75 | acc_intro < 68 | annot == "GeneMark-ES-F" | annot == "VDAG_Jr2_Annotation" | annot == "CodingQuarry", as.character(int_gene$annot),"")),
    fontface = 'plain',  color ="black", size= 2,
    box.padding= unit(0.35, "lines"),
    point.padding= unit(0.5, "lines"),
    force = 2, 
    segment.color= 'grey50') +
  geom_smooth(method=lm) +
  theme_bw() +
  theme(legend.position = 'none')


ggsave(file="intron_gene_corr.pdf", int_gene2, width = 120, height = 120, units = "mm")



