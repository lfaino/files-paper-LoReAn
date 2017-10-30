setwd("C:/Users/luigi/Dropbox/LoReAn_paper/Data/files-paper-LoReAn/V_dahliae/newGtf/orthomcl")

library(ggplot2)
library(gridExtra)
library(scales)
library(agricolae)

allData <- read.table("all.ready.txt", header=F)
colnames(allData) <- c("Pipeline", "Cov", "Sin_Len")


boxCov2 <- ggplot(allData, aes(x=Pipeline, y=Cov, fill= Pipeline)) + 
  geom_boxplot( lwd = 0.6, outlier.shape = NA) +
  #geom_dotplot(binaxis= "y", stackdir= "center", dotsize=0.35, lwd = 0.3) +
  geom_jitter(shape = 20, position=position_jitter(0.2), size = 0.5) +
  scale_fill_manual(values=c("#658361", "#2C9E4B", "#9ED763", "#FFF2B2")) +
  theme_classic() + labs(y = "SingletonCoverage") +
  theme(axis.title.x = element_blank(), legend.position='none')

denCov <- ggplot(allData, aes(x=Cov, fill= Pipeline)) + 
  geom_density(alpha = 0.3) +
  scale_fill_manual(values=c("#658361", "#2C9E4B", "#9ED763", "#FFF2B2"))

#Using faceting, interesting but not best
denCov2 <- ggplot(allData, aes(x=Cov)) + 
  geom_density(alpha = 0.3) +
  scale_fill_manual(values=c("#658361", "#2C9E4B", "#9ED763", "#FFF2B2")) +
  facet_grid(Pipeline ~.)


boxLen2 <- ggplot(allData, aes(x=Pipeline, y=log2(Sin_Len), fill= Pipeline)) + 
  geom_boxplot( lwd = 0.6, outlier.shape = NA) +
  #geom_dotplot(binaxis= "y", stackdir= "center", dotsize=0.35, lwd = 0.3) +
  geom_jitter(shape = 20, position=position_jitter(0.2), size = 0.5) +
  scale_fill_manual(values=c("#658361", "#2C9E4B", "#9ED763", "#FFF2B2")) +
  theme_classic() + labs(y = "SingletonLength") +
  theme(axis.title.x = element_blank(), legend.position='none') 
  #ylim(c(0,5000))

denLen <- ggplot(allData, aes(x=Sin_Len, fill= Pipeline)) + 
  geom_density(alpha = 0.3) +
  scale_fill_manual(values=c("#658361", "#2C9E4B", "#9ED763", "#FFF2B2")) +
  xlim(c(0,5000))

both <- grid.arrange(boxCov2, boxLen2, ncol=2, nrow=1)
ggsave(file=paste0("./singleton_ColLenOrthp.2.pdf"), both, width = 90, height = 50, units = "mm")
dev.off()

av <- aov(Cov ~ Pipeline, data = allData)
HSD.test(av,"Pipeline", console=TRUE)

library(doBy)
Covdata <- summaryBy(Cov ~ Pipeline, data = allData, FUN=c(median, mean, sd))
Lendata <- summaryBy(Sin_Len ~ Pipeline, data = allData, FUN=c(median, mean, sd))





