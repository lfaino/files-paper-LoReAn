#setwd("~/Dropbox/LoReAn_paper/Data/ReDo_JR2")
setwd("C:/Users/luigi/Dropbox/LoReAn_paper/Data/files-paper-LoReAn/V_dahliae/newGtf/plots/")

library(reshape2)
library(ggplot2)
library(gridExtra)
library(agricolae)


allData <- read.csv("../SummaryExonerate1.csv", header=T)
exclude <- c("maker", "coding", "augustus" , "braker", "genemark-es", "genemark-es-fu")
data48 <- subset(allData, !(Pipeline %in% exclude))
data48 <- droplevels(data48)

ex_gene48 <- data48  [ , c(1, 6:8, 40:42)]
ov_gene48 <- data48  [ , c(1, 11:13, 40:42)]
ex_tran48 <- data48  [ , c(1, 19:21, 40:42)]
ov_tran48 <- data48  [ , c(1, 24:26, 40:42)]
ex_exon48 <- data48  [ , c(1, 32:34, 40:42)]
ov_exon48 <- data48  [ , c(1, 37:42)]

melt_ex_gene48  <- melt(ex_gene48 , id.vars=c("Annotation", "masking", "abOption", "Pipeline"))
melt_ov_gene48  <- melt(ov_gene48 , id.vars=c("Annotation", "masking", "abOption", "Pipeline"))
melt_ex_tran48  <- melt(ex_tran48 , id.vars=c("Annotation", "masking", "abOption", "Pipeline"))
melt_ov_tran48  <- melt(ov_tran48 , id.vars=c("Annotation", "masking", "abOption", "Pipeline"))
melt_ex_exon48  <- melt(ex_exon48 , id.vars=c("Annotation", "masking", "abOption", "Pipeline"))
melt_ov_exon48  <- melt(ov_exon48 , id.vars=c("Annotation", "masking", "abOption", "Pipeline"))

df <- c("melt_ex_gene48", "melt_ov_gene48", "melt_ex_tran48", "melt_ov_tran48", "melt_ex_exon48", "melt_ov_exon48")
var = c("masking", "abOption", "Pipeline")
mydfs <- lapply(df, get) #This makes an actual list of the dataframes
names(mydfs) <- df

for (i in df){
  mydfs[[i]]$masking <- relevel(mydfs[[i]]$masking, "partMasked")
  mydfs[[i]]$masking <- relevel(mydfs[[i]]$masking, "Masked")
  mydfs[[i]]$abOption <- relevel(mydfs[[i]]$abOption, "Fungus")
  mydfs[[i]]$abOption <- relevel(mydfs[[i]]$abOption, "Braker")
}

#exact <- c("melt_ex_gene48", "melt_ex_tran48", "melt_ex_exon48")
#over <- c("melt_ov_gene48", "melt_ov_tran48", "melt_ov_exon48")
#if (i %in% exact){
#  a = c(94, 96, 98, 100)
#} else {
#  a = limits=c(35,45,55,65)
#}

for (i in df){
  mask <- ggplot(mydfs[[i]], aes(x=variable, y=value, fill= masking)) + 
    geom_boxplot( lwd = 0.3) +
    geom_dotplot(binaxis= "y", stackdir= "center", position=position_dodge(0.8), dotsize=0.35, lwd = 0.3) +
    scale_fill_manual(values=c("#339999", "#FFFFCC", "#99CCCC")) +
    theme_classic() + labs(y = i) +
    theme(axis.title.x = element_blank(), legend.position='none')
  
  mask2 <- mask + scale_y_continuous(breaks=pretty(mydfs[[i]]$value, n=5))
   
  ab <- ggplot(mydfs[[i]], aes(x=variable, y=value, fill=abOption)) + 
    geom_boxplot( lwd = 0.3) +
    geom_dotplot(binaxis= "y", stackdir= "center", position=position_dodge(0.8), dotsize=0.35, lwd = 0.3) +
    scale_fill_manual(values=c("#658361", "#2C9E4B", "#9ED763", "#FFF2B2")) +
    theme_classic() + labs(y = i) +
    theme(axis.title.x = element_blank(), legend.position="none")
  ab2 <- ab + scale_y_continuous(breaks=pretty(mydfs[[i]]$value, n=5))
  
  pplot <- ggplot(mydfs[[i]], aes(x=variable, y=value, fill= Pipeline)) +
    geom_boxplot( lwd = 0.3) + 
    geom_dotplot(binaxis= "y", stackdir= "center", position=position_dodge(0.8), dotsize=0.35, lwd = 0.3) +
    scale_fill_manual(values=c("#FAAF08", "#FA812F", "#FA4032", "#FEF3E2")) +
    theme_classic() + labs(y = i) +
    theme(axis.title.x = element_blank(), legend.position="none")
  pplot2 <- pplot + scale_y_continuous(breaks=pretty(mydfs[[i]]$value, n=5))
  
  test <- grid.arrange(mask2, ab2, pplot2, ncol=3, nrow=1)
  ggsave(file=paste0("./",i,"QualSummaryPlots.pdf"), test, width = 170, height = 80, units = "mm")
  dev.off()
}

melt_ex_geneAll <- melt(allData [ , c(1, 6:8)], id.vars = "Annotation")


melt_ex_geneAll$Annotation <- factor(
  melt_ex_geneAll$Annotation, levels=c("augustus", "Braker", "GeneMArks_Fu", "GeneMArks", "CodingQuarry", "MAKER2", "BAP_MGS_Br", "BAPplus_MGS_Br", "LoReAn_MGnS_Br", "LoReAn_MGS_Br", "BAP_MGS_Br_Fu", "BAPplus_MGS_Br_Fu", "LoReAn_MGnS_Br_Fu", "LoReAn_MGS_Br_Fu", "BAP_MGS_Fu", "BAPplus_MGS_Fu", "LoReAn_MGnS_Fu", "LoReAn_MGS_Fu", "BAP_MGS", "BAPplus_MGS", "LoReAn_MGnS", "LoReAn_MGS", "BAP_FGS_Br", "BAPplus_FGS_Br", "LoReAn_FGnS_Br", "LoReAn_FGS_Br", "BAP_FGS_Br_Fu", "BAPplus_FGS_Br_Fu", "LoReAn_FGnS_Br_Fu", "LoReAn_FGS_Br_Fu", "BAP_FGS_Fu", "BAPplus_FGS_Fu", "LoReAn_FGnS_Fu", "LoReAn_FGS_Fu", "BAP_FGS", "BAPplus_FGS", "LoReAn_FGnS", "LoReAn_FGS", "BAP_pMGS_Br", "BAPplus_pMGS_Br", "LoReAn_pMGnS_Br", "LoReAn_pMGS_Br", "BAP_pMGS_Br_Fu", "BAPplus_pMGS_Br_Fu", "LoReAn_pMGnS_Br_Fu", "LoReAn_pMGS_Br_Fu", "BAP_pMGS_Fu", "BAPplus_pMGS_Fu", "LoReAn_pMGnS_Fu", "LoReAn_pMGS_Fu", "BAP_pMGS", "BAPplus_pMGS", "LoReAn_pMGnS", "LoReAn_pMGS"))
            
            
            
melt_ex_geneAll$value2 <- melt_ex_geneAll$value*100

ex_gene_allD <- ggplot(melt_ex_geneAll, aes(x=value2, y=Annotation, colour=variable)) +
  geom_point(size = 2) +
  #This is used to make the horizontal lines at each y-value
  geom_segment(aes(x=25, y=as.numeric(Annotation), xend=75, yend=as.numeric(Annotation)), colour="grey") +
  scale_x_discrete(limits = c(25, 35, 45, 55, 65, 75)) +
  theme_classic() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position="none") 
ggsave(file=("./ex_gene_allD_horizontal.pdf"), width = 85, height = 160, units = "mm")

#Comparisons, summary, stats


ex_gene48_2 <- data48  [ , c(1, 2, 6:8, 40:42)]
ov_gene48_2 <- data48  [ , c(1, 2, 11:13, 40:42)]
ov_tran48_2 <- data48  [ , c(1, 14, 24:26, 40:42)]
ov_exon48_2 <- data48  [ , c(1, 27, 37:42)]

ex_tran48_2 <- data48  [ , c(1, 14, 16, 19:21, 40:42)]
ex_exon48_2 <- data48  [ , c(1, 27, 29, 32:34, 40:42)]


ex_gene48_3 <- data48  [ , c(1, 2, 4, 6:8, 40:42)]
ex_tran48_3 <- data48  [ , c(1, 14, 17, 19:21, 40:42)]
ex_exon48_3 <- data48  [ , c(1, 27, 30, 32:34, 40:42)]

df2 <- c( "ex_gene48_2", "ov_gene48_2", "ov_tran48_2", "ov_exon48_2")
mydf2 <- lapply(df2, get)
names(mydf2) <- df2
for (i in df2){
  colnames(mydf2[[i]]) <- c("Annotation", "PredCount", "Specificity",
                            "Sensitivity", "Accuracy", "masking", 
                            "abOption", "Pipeline")
  }

df3 <- c( "ex_tran48_2", "ex_exon48_2")
mydf3 <- lapply(df3, get)
names(mydf3) <- df3
for (i in df3){
  colnames(mydf3[[i]]) <- c("Annotation", "PredCount", "AvgLength", "Specificity",
                            "Sensitivity", "Accuracy", "masking", 
                            "abOption", "Pipeline")
  }

df4 <- c("ex_gene48_3", "ex_tran48_3", "ex_exon48_3")
mydf4 <- lapply(df4, get)
names(mydf4) <- df4
for (i in df4){
  colnames(mydf4[[i]]) <- c("Annotation", "PredCount", "ExactMatchRef",
                            "Specificity", "Sensitivity", "Accuracy", "masking", 
                            "abOption", "Pipeline")
}



var2 = c("PredCount", "Specificity", "Sensitivity", "Accuracy") #This is setting the columns to compare Prediction, Spec, Sens, Acc
var3 = c("PredCount", "AvgLength", "Specificity", "Sensitivity", "Accuracy")
var4 = c("PredCount", "ExactMatchRef", "Specificity", "Sensitivity", "Accuracy")

compare = c("masking", "abOption", "Pipeline") # This is to go through HSD comparing the various options

#This can be run for df2 and df3, but the var list is different. Just change the df# to get 2 or 3.
for (i in df3){ # for each of column Spec, Sens, Acc
  for (z in var3){
    av <- aov(mydf3[[i]][ ,z] ~ masking + abOption + Pipeline, data = mydf3[[i]])
    sink(paste0("./LoReAn.stats.",i,".txt"), append=T)
    print(paste("Start for",colnames(mydf3[[i]][z]),"Variable"))
    for (q in compare){ # for each of Masking, abOption, pipe
      print("*************************")
      (HSD.test(av, q, console=TRUE)) # 1st run = HSD.test(av, Masking, console=TRUE)
      print("*************************")
    } 
    print ("**************************")
  }
  sink()
}
closeAllConnections()


for (i in df4){ # for each of column Spec, Sens, Acc
  for (z in var4){
    av <- aov(mydf4[[i]][ ,z] ~ masking + abOption + Pipeline, data = mydf4[[i]])
    sink(paste0("./LoReAn.stats.2.",i,".txt"), append=T)
    print(paste("Start for",colnames(mydf4[[i]][z]),"Variable"))
    for (q in compare){ # for each of Masking, abOption, pipe
      print("*************************")
      (HSD.test(av, q, console=TRUE)) # 1st run = HSD.test(av, Masking, console=TRUE)
      print("*************************")
    } 
    print ("**************************")
  }
  sink()
}
closeAllConnections()

