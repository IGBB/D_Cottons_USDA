
library("ggplot2")
library("dplyr")
library("reshape2")
library("gridExtra")
library("data.table")

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


Species_translate <- as.data.frame(c(levels(DF$qspecies), "D9-4"))
colnames(Species_translate) <- c("Species")
Species_translate$Species <- as.character(Species_translate$Species)
Species_translate$latin <- c("G. thurberi ", "G. turneri", "G. schwendemanii","G. armourianum", "G. harknessii", "G. davidsonii","G. klotzschianum", "G. aridum", "G. raimondii","G. gossypioides", "G. lobatum", "G. trilobum", "G. laxum")
Species_translate$order <- c("I", "D", "H", "B", "C", "K", "L", "E", "A", "M", "G", "J", "F")

######
#Table
######

groups <- DF %>% group_by(`sspecies`, `qspecies`)

#dS Table
dS.summary.table <- summarize(groups, LQ = quantile(dS, 0.25), Median = median(dS), UQ = quantile(dS, 0.75))
levels(dS.summary.table$sspecies) <- c(levels(dS.summary.table$sspecies), "D1-35", levels(Species_translate$order))
levels(dS.summary.table$qspecies) <- c(levels(dS.summary.table$qspecies), "D9-4", levels(Species_translate$order))
dS.summary.table[nrow(dS.summary.table) + 1,] <- list("D1-35", "D1-35", NA, NA, NA)
dS.summary.table[nrow(dS.summary.table) + 1,] <- list("D9-4", "D9-4", NA, NA, NA)
dS.summary.table$sspecies <- Species_translate$order[match(dS.summary.table$sspecies, Species_translate$Species)]
dS.summary.table$qspecies <- Species_translate$order[match(dS.summary.table$qspecies, Species_translate$Species)]
dS.summary.table$IQR <- paste(format(dS.summary.table$LQ, digits = 3), format(dS.summary.table$UQ, digits = 3), sep="-")
dS.summary.table$IQR <- paste("(", dS.summary.table$IQR,")", sep="")
dS.summary.table$concat <- paste(format(dS.summary.table$Median, digits = 3), dS.summary.table$IQR, sep = " ")


dS.summary.df <- as.matrix(t(acast(dS.summary.table, sspecies~qspecies, value.var = "concat")))
dS.summary.df <- dS.summary.df[order(rownames(dS.summary.df)),order(colnames(dS.summary.df))]
dS.summary.df <- as.data.frame(replace(dS.summary.df, is.na(dS.summary.df), t(dS.summary.df)[is.na(dS.summary.df)]))
dS.summary.df[lower.tri(dS.summary.df, diag = T)] <- NA

#dN Table

dN.summary.table <- summarize(groups, LQ = quantile(dN, 0.25), Median = median(dN), UQ = quantile(dN, 0.75))
levels(dN.summary.table$sspecies) <- c(levels(dN.summary.table$sspecies), "D1-35", levels(Species_translate$order))
levels(dN.summary.table$qspecies) <- c(levels(dN.summary.table$qspecies), "D9-4", levels(Species_translate$order))
dN.summary.table[nrow(dN.summary.table) + 1,] <- list("D1-35", "D1-35", NA, NA, NA)
dN.summary.table[nrow(dN.summary.table) + 1,] <- list("D9-4", "D9-4", NA, NA, NA)
dN.summary.table$sspecies <- Species_translate$order[match(dN.summary.table$sspecies, Species_translate$Species)]
dN.summary.table$qspecies <- Species_translate$order[match(dN.summary.table$qspecies, Species_translate$Species)]
dN.summary.table$IQR <- paste(format(dN.summary.table$LQ, digits = 2), format(dN.summary.table$UQ, digits = 3), sep="-")
dN.summary.table$IQR <- paste("(", dN.summary.table$IQR, ")", sep="")
dN.summary.table$concat <- paste(format(dN.summary.table$Median, digits = 2), dN.summary.table$IQR, sep = " ")


dN.summary.df <- as.matrix(acast(dN.summary.table, sspecies~qspecies, value.var = "concat"))
dN.summary.df <- dN.summary.df[order(rownames(dN.summary.df)),order(colnames(dN.summary.df))]
dN.summary.df <- as.data.frame(replace(dN.summary.df, is.na(dN.summary.df), t(dN.summary.df)[is.na(dN.summary.df)]))
dN.summary.df[upper.tri(dN.summary.df, diag = T)] <- NA



#######
#Combine tables
#######

dN.summary.df <- as.data.frame(lapply(dN.summary.df, as.character), stringsAsFactors=FALSE)
dS.summary.df <- as.data.frame(lapply(dS.summary.df, as.character), stringsAsFactors=FALSE)

Total.df <- replace(dS.summary.df, is.na(dS.summary.df), dN.summary.df[is.na(dS.summary.df)])

colnames(Total.df) <- Species_translate$latin[match(colnames(Total.df), Species_translate$order)]
row.names(Total.df) <- colnames(Total.df)




write.table(Total.df, "dNdS_table.txt", quote=F, sep="\t")

