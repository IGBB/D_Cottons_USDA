
library("ggplot2")
library("dplyr")
library("reshape2")

######
#Read Data
######
setwd("~/GitHub_Repos/D_Cottons_USDA/analysis/dNdS_CNV")
DF <- read.csv("final_dNdS.csv", header = TRUE)



######
#Pruning
######


DF$`Species Comparison` <- as.factor(paste(DF$qspecies, DF$sspecies, sep = " vs. "))

filtered <- !(DF$X %in% DF[DF$dS > 0.6,]$X)
DF <- DF[filtered,]




######
#Table
######

groups <- DF %>% group_by(`sspecies`, `qspecies`)
dS.summary.table <- summarize(groups, LQ = quantile(dS, 0.25), Median = median(dS), UQ = quantile(dS, 0.75))
levels(dS.summary.table$sspecies) <- c(levels(dS.summary.table$sspecies), "D1-35")
levels(dS.summary.table$qspecies) <- c(levels(dS.summary.table$qspecies), "D9-4")
dS.summary.table[nrow(dS.summary.table) + 1,] <- c("D1-35", "D1-35", NA, NA, NA)
dS.summary.table[nrow(dS.summary.table) + 1,] <- c("D9-4", "D9-4", NA, NA, NA)
dS.summary.table$IQR <- paste(dS.summary.table$LQ, dS.summary.table$UQ, sep="-")
dS.summary.table$IQR <- paste("(", dS.summary.table$IQR, ")", sep="")
dS.summary.table$concat <- paste(dS.summary.table$Median, dS.summary.table$IQR, sep = " ")
dS.summary.df <- as.data.frame(t(acast(dS.summary.table, sspecies~qspecies, value.var = "concat")))
dS.summary.df <- dS.summary.df[order(rownames(dS.summary.df)),order(colnames(dS.summary.df))]
dS.summary.df[dS.summary.df=="NA (NA-NA)"] <- NA


dN.summary.table <- summarize(groups, LQ = quantile(dN, 0.25), Median = median(dN), UQ = quantile(dN, 0.75))
levels(dN.summary.table$sspecies) <- c(levels(dN.summary.table$sspecies), "D1-35","---")
levels(dN.summary.table$qspecies) <- c(levels(dN.summary.table$qspecies), "D9-4","---")
dN.summary.table[nrow(dN.summary.table) + 1,] <- c("D1-35", "D1-35", NA, NA, NA)
dN.summary.table[nrow(dN.summary.table) + 1,] <- c("D9-4", "D9-4", NA, NA, NA)
dN.summary.table$IQR <- paste(dN.summary.table$LQ, dN.summary.table$UQ, sep="-")
dN.summary.table$IQR <- paste("(", dN.summary.table$IQR, ")", sep="")
dN.summary.table$concat <- paste(dN.summary.table$Median, dN.summary.table$IQR, sep = " ")
dN.summary.df <- as.data.frame(acast(dN.summary.table, sspecies~qspecies, value.var = "concat"))
dN.summary.df <- dN.summary.df[order(rownames(dN.summary.df)),order(colnames(dN.summary.df))]
dN.summary.df[dN.summary.df=="NA (NA-NA)"] <- NA

dN.summary.df <- data.frame(lapply(dN.summary.df, as.character), stringsAsFactors=FALSE)
dS.summary.df <- data.frame(lapply(dS.summary.df, as.character), stringsAsFactors=FALSE)

dN.summary.df[is.na(dN.summary.df)] <- "---" #THIS DOESN"T FREAKING WORK AND I DON"T KNOW WHY 
Total.df <- replace(dS.summary.df, is.na(dS.summary.df), dN.summary.df[is.na(dS.summary.df)])


Df.melted <- melt(DF, id.vars = "Species Comparison", measure.vars = c("dS", "dN"))
colnames(Df.melted)[3] <- c("Synonymous Substitution Rate")

newDf.melted <- Df.melted[Df.melted$variable == "dS",]

######
#Plotting
######

#png("Figure_dNdSboxplot.png", 5000, 5000, pointsize=12, res=600)
#ggplot(Df.melted, aes(x=`Species Comparison`, y = `Synonymous Substitution Rate`, fill=variable)) + geom_boxplot(position = position_dodge(1), outlier.shape = NA) + 
#  theme(legend.position = c(0.95,0.9), legend.title = element_blank(), axis.title.x = element_blank(), axis.text.x = element_text(color = "black", size = 13)) + 
#  scale_y_continuous(limits = c(0,0.175)) +
#  stat_summary(fun.y=median, geom="point", size=1, color="white", position = position_dodge(1)) + scale_fill_manual(values = c("red3", "black"))
#dev.off()




png("Figure_dSDistribution.png", 5000, 5000, pointsize=12, res=600)
ggplot(newDf.melted, aes(newDf.melted$`Synonymous Substitution Rate`, colour = newDf.melted$`Species Comparison`)) + 
  geom_density(adjust = 1/5) + 
  theme(legend.position = c(0.6, 0.6), axis.title.y = element_blank(), legend.title = element_blank()) + 
  xlab("Synonymous Substitution Rate") 
dev.off()




#png("Figure_Malvaceae_Divergence.png", 5000, 5000, pointsize=12, res=600)
#ggplot(Malvaceae_DF, aes(Malvaceae_DF$dS)) + 
#  geom_density(adjust = 1/5) + geom_vline(xintercept = median(Malvaceae_DF$dS)) +
#  theme(axis.title.y = element_blank(), legend.title = element_blank()) + 
#  xlab("Synonymous Substitution Rate") 
#dev.off()

