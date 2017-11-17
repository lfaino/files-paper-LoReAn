setwd("C:/Users/luigi/Dropbox/LoReAn_paper/Data/ReDo_JR2/")
library(reshape2)
library(ggplot2)
library(gridExtra)
library(agricolae)
library(ggrepel)

##############################################
######## Introns #############################
JR2_intron <- read.csv("./../files-paper-LoReAn/V_dahliae/intron//onlyIntrons.csv", header=T, fill=T)


acclines2 <- data.frame(
  x = c(140, 160, 180),
  y = c(0, 0, 0),
  xend = c(0, 0, 0),
  yend = c(140, 160, 180),
  color = c("black", "red", "orange")
) 

intron.jr <- ggplot(JR2_intron, aes(x=spe_intro, y=sens_intro)) +
  geom_point(aes(color = annot), size =1.25) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend, color=color), linetype="dashed", data=acclines2) +
  coord_cartesian(ylim=c(60, 85), xlim=c(60, 80)) +
  geom_label_repel(
    aes(fill = JR2_intron$annot, label= ifelse(acc_intro > 75 | acc_intro < 68 | annot == "GeneMArks_Fu.gtf.introns" | annot == "VDAG_Jr2_Annotation.v5.introns" | annot == "CodingQuarry_out_PredictedPass.gff3.introns", as.character(JR2_intron$annot),"")),
    fontface = 'plain',  color ="white", size= 2.5,
    box.padding= unit(0.35, "lines"),
    point.padding= unit(0.5, "lines"),
    force = 2, 
    segment.color= 'grey50') +
  theme_bw() +
  theme(legend.position = 'none')

ggsave(file="jr2.introns.pdf", intron.jr, width = 120, height = 120, units = "mm")



###################################
### Intron Analysis  ##############


int_gene <- read.csv2("../files-paper-LoReAn/V_dahliae/newGtf/Summary_1_gene_intron.csv", header=T,
                      sep= ",", dec = ".", stringsAsFactors = F)
int_gene$value <- paste (int_gene$genome, "_", int_gene$pipeline)
int_gene$accuracy1 <- int_gene$accuracy * 100

#Scatter plot with regression lines
int_gene_plot <- ggplot(int_gene, aes(x=acc_intro, y=accuracy1)) +
  geom_point() +
  xlab("Intron Accuracy") +
  ylab("Exact Gene Accuracy") +
  geom_smooth(method=lm) +
  theme_bw()

with(int_gene, cor(as.numeric(acc_intro), as.numeric(accuracy1)), method="spearman")
with(int_gene, cor.test(as.numeric(acc_intro), as.numeric(accuracy1), method=c("spearman")))

plot(int_gene$acc_intro,int_gene$accuracy1,
     xlab="Intron Accuracy", ylab="Exact Gene Accuracy")
abline(lm(int_gene$accuracy1 ~ int_gene$acc_intro))

int_gene2 <- ggplot(int_gene, aes(x=acc_intro, y=accuracy1)) +
  geom_point(aes(color = pipeline), size =1.25) +
  #geom_segment(aes(x=x, y=y, xend=xend, yend=yend, color=color), linetype="dashed", data=acclines_r) +
  #coord_cartesian(ylim=c(50, 100), xlim=c(50, 100)) +
  geom_label_repel(
    aes(fill = pipeline, label= annot),
    fontface = 'plain',  color ="white", size= 1.5,
    box.padding= unit(0.35, "lines"),
    point.padding= unit(0.5, "lines"),
    force = 2, 
    segment.color= 'grey50') +
  geom_smooth(method=lm) +
  theme_bw() +
  theme(legend.position = 'none')


ggsave(file="intron_gene_corr.pdf", int_gene2, width = 120, height = 120, units = "mm")


