setwd("C:/Users/luigi/Dropbox/LoReAn_paper/Data/files-paper-LoReAn/V_dahliae/newGtf/orthomcl")
library(ggplot2)
library(gridExtra)
library(scales)
library(ggpubr)

allData1 <- read.csv("AllPie.cov.exonsPIE.csv", header=T)
row.names(allData1) <- allData1$X
allData1$X <- NULL
allData <- as.data.frame(t(allData1))

pdf(file="maker.pdf")

ggplot(allData, aes(x="", y = allData$maker, fill = row.names(allData))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), axis.text.x= element_blank()) +
  geom_text(aes(y = (allData$maker / 2) + c(0, cumsum(allData$maker)[-length(allData$maker)]),
                label = percent(allData$maker/ sum(allData$maker))), size = 10, color = "#FFFFFF")
dev.off()
pdf(file="lorean.pdf")

ggplot(allData, aes(x="", y = allData$lorean, fill = row.names(allData))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), axis.text.x= element_blank()) +
  geom_text(aes(y = (allData$lorean / 2) + c(0, cumsum(allData$lorean)[-length(allData$lorean)]),
                label = percent(allData$lorean / sum(allData$lorean))), size = 10, color = "#FFFFFF")
dev.off()
pdf(file="bap.pdf")

ggplot(allData, aes(x="", y = allData$bap, fill = row.names(allData))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), axis.text.x= element_blank()) +
  geom_text(aes(y = (allData$bap / 2) + c(0, cumsum(allData$bap)[-length(allData$bap)]),
                label = percent(allData$bap / sum(allData$bap))), size = 10, color = "#FFFFFF")
dev.off()

pdf(file="cq.pdf")
ggplot(allData, aes(x="", y = allData$cq, fill = row.names(allData))) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  theme_minimal() +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        panel.grid = element_blank(), axis.text.x= element_blank()) +
  geom_text(aes(y = (allData$cq / 2) + c(0, cumsum(allData$cq)[-length(allData$cq)]),
                label = percent(allData$cq / sum(allData$cq))), size = 10, color = "#FFFFFF")
dev.off()
