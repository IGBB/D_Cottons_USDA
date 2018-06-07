library(reshape2)

df <- read.table("D.kmer.table", sep="\t", stringsAsFactors=F)
df$V3 <- gsub(" bp|%|,", "", df$V3)
df$V1 <- gsub("D3_", "D3", df$V1)

cast <- as.data.frame(acast(df, V1~V2))

GS <- read.table("genome_sizes_combine.txt", header=T, sep="\t", stringsAsFactors=F)
GS$Species <- gsub("\\.|-", "_", GS$Species)

GSrm <- setdiff(GS$Species, row.names(cast))

row.names(GS) <- GS$Species

GS <- subset(GS, !row.names(GS) %in% GSrm)

tog <- merge(cast, GS, by="row.names")
names(tog) <- gsub(" ", "_", names(tog))

tog$HapSize <- round(as.numeric(as.character(tog$Genome_Haploid_Length))/1000000,0)
row.names(tog) <- tog$Row.names

png("GSvsKmerGS.png", 5000, 5000, pointsize=12, res=600)
ggplot(tog,aes(Size,HapSize)) + geom_smooth(method='lm')+geom_point(aes(color=Species)) + ylim(500,1000) + xlim(500,1000)+
     stat_poly_eq(formula = tog$HapSize ~ tog$Size, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)

dev.off()

TE <- read.table("TE.totals.tbl", header=T, row.names=1, sep="\t", stringsAsFactors=F)
row.names(TE) <- gsub("\\.|-", "_", row.names(TE))

TErm <- setdiff(row.names(TE), row.names(tog))
TE <- subset(TE, !row.names(TE) %in% TErm)

togTE <- merge(tog, TE, by="row.names")
row.names(togTE) <- togTE$Row.names

togTE$Row.names <- NULL
togTE$Row.names <- NULL

togTE$TEtotal <- round(togTE$TEtotal, 0)
togTE$RepSize <- round(as.numeric(as.character(tog$Genome_Repeat_Length))/1000000,0)

png("GSvsKmerTE.png", 5000, 5000, pointsize=12, res=600)
ggplot(togTE,aes(TEtotal,RepSize)) + geom_smooth(method='lm')+geom_point(aes(color=Species))  +
   stat_poly_eq(formula = togTE$RepSize ~ togTE$TEtotal, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
dev.off()



png("ExpTEvsKmerGS.png", 5000, 5000, pointsize=12, res=600)
ggplot(togTE,aes(HapSize,TEtotal)) + geom_smooth(method='lm')+geom_point(aes(color=Species))  +
   stat_poly_eq(formula = togTE$TEtotal ~ togTE$HapSize, aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), parse = TRUE)
dev.off()



