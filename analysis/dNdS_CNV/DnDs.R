
library("ggplot2")
library("dplyr")
library("reshape2")
library("gridExtra")

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
dS.summary.table$concat <- paste(dS.summary.table$Median, dS.summary.table$IQR, sep = "\t")
dS.summary.df <- as.data.frame(t(acast(dS.summary.table, sspecies~qspecies, value.var = "concat")))
dS.summary.df <- dS.summary.df[order(rownames(dS.summary.df)),order(colnames(dS.summary.df))]
dS.summary.df[dS.summary.df=="NA\t(NA-NA)"] <- NA


dN.summary.table <- summarize(groups, LQ = quantile(dN, 0.25), Median = median(dN), UQ = quantile(dN, 0.75))
levels(dN.summary.table$sspecies) <- c(levels(dN.summary.table$sspecies), "D1-35","---")
levels(dN.summary.table$qspecies) <- c(levels(dN.summary.table$qspecies), "D9-4","---")
dN.summary.table[nrow(dN.summary.table) + 1,] <- c("D1-35", "D1-35", NA, NA, NA)
dN.summary.table[nrow(dN.summary.table) + 1,] <- c("D9-4", "D9-4", NA, NA, NA)
dN.summary.table$IQR <- paste(dN.summary.table$LQ, dN.summary.table$UQ, sep="-")
dN.summary.table$IQR <- paste("(", dN.summary.table$IQR, ")", sep="")
dN.summary.table$concat <- paste(dN.summary.table$Median, dN.summary.table$IQR, sep = "\t")
dN.summary.df <- as.data.frame(acast(dN.summary.table, sspecies~qspecies, value.var = "concat"))
dN.summary.df <- dN.summary.df[order(rownames(dN.summary.df)),order(colnames(dN.summary.df))]
dN.summary.df[dN.summary.df=="NA\t(NA-NA)"] <- NA

dN.summary.df <- data.frame(lapply(dN.summary.df, as.character), stringsAsFactors=FALSE)
dS.summary.df <- data.frame(lapply(dS.summary.df, as.character), stringsAsFactors=FALSE)

dN.summary.df[is.na(dN.summary.df)] <- "---" #THIS DOESN"T FREAKING WORK AND I DON"T KNOW WHY 
Total.df <- replace(dS.summary.df, is.na(dS.summary.df), dN.summary.df[is.na(dS.summary.df)])
colnames(Total.df) <- c("G. thurberi ", "G. turneri", "G. schwendemanii","G. armourianum", "G. harknessii", "G. davidsonii","G. klotzschianum", "G. aridum", "G. raimondii","G. gossypioides", "G. lobatum", "G. trilobum", "G. laxum")
row.names(Total.df) <- colnames(Total.df)

Df.melted <- melt(DF, id.vars = "Species Comparison", measure.vars = c("dS", "dN"))
colnames(Df.melted)[3] <- c("Synonymous Substitution Rate")

newDf.melted <- Df.melted[Df.melted$variable == "dS",]

######
#Plotting
######


#png("Table_dSDistribution.png", 5000, 5000, pointsize=12, res=600)
#ggplot(newDf.melted, aes(newDf.melted$`Synonymous Substitution Rate`, colour = newDf.melted$`Species Comparison`)) + 
#  geom_density(adjust = 1/5) + 
#  theme(legend.position = c(0.6, 0.6), axis.title.y = element_blank(), legend.title = element_blank()) + 
#  xlab("Synonymous Substitution Rate") 
#dev.off()


write.table(Total.df, "dNdS_table.txt")

