library(ggplot2)
library(scales)
library(ggrepel)
library(factoextra)
library(reshape2)
library(gridExtra)
library(geomorph)
library(testit)
library(ape)
library(geiger)
library(phytools)

#############################################
sessionInfo()
# R version 3.3.3 (2017-03-06)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 14393)
# 
# locale:
# [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
# [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] reshape2_1.4.2       gridExtra_2.2.1      factoextra_1.0.4     ggrepel_0.6.5        plotrix_3.6-4        BiocInstaller_1.22.3 scales_0.4.1         ggplot2_2.2.1       
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.10     digest_0.6.12    ggpubr_0.1.2     grid_3.3.3       plyr_1.8.4       gtable_0.2.0     magrittr_1.5     stringi_1.1.5    lazyeval_0.2.0   labeling_0.3     tools_3.3.3     
# [12] stringr_1.2.0    munsell_0.4.3    colorspace_1.3-2 tibble_1.3.0    
# 
##############################################

# import file with all counts
# evaluate if 0.01% cutoff is reasonable (here, cluster 274)

data <- read.table("comparative_analysis_counts.annotated.txt", header = T, row.names=1, sep="\t")
data$size <- as.numeric(rowSums(data[,-1]))
data$percent <- cumsum(data$size)/sum(data$size)

png("cotton.cutoff.png", 5000, 5000, pointsize=12, res=600)
ggplot(data, aes(x=cluster, y=percent)) + geom_line(size=1) + geom_vline(xintercept=274, color='yellow3', size=1) + scale_x_log10(labels=comma) + scale_y_log10() + geom_vline(xintercept=0, color="grey")+ geom_hline(yintercept=0, color="grey")
dev.off()


################### ordination ###################
### PCoA of scaled data ###

annot_clust <- read.table("annotated.counts.txt", header = T, row.names=1, sep="\t")

ord_table <- annot_clust[,(3:34)]
Mb_table <- ord_table*0.0095 
# 0.0095 multiplier represents # reads (x) * 95nt/read * 1 kb/1000nt * 1Mb/1000kb * 100% = # reads * 0.0095 = # Mb in entire genome for that cluster 

# convert to percent of genome to balance the numbers for ordination, asks are the relative proportions the same
#genome	size	reads
#D01	841	98942
#D02.1	856	100706
#D02.2	910	107059
#D10	910	107059
#D11	929	109295
#D3D	910	107059
#D3K	880	103530
#D4	919	108118
#D5	880	103530
#D6	841	98942
#D7	934	109883
#D8	851	100118
#D9	934	109883
#A1	1667	196118
#A2	1698	199765
#F1	1311	154236
perc_table <- data.frame( lapply(Mb_table[1:2], function(x) x/841), 
lapply(Mb_table[3], function(x) x/856), 
lapply(Mb_table[4], function(x) x/910),
lapply(Mb_table[24:26], function(x) x/910),
lapply(Mb_table[27], function(x) x/929),
lapply(Mb_table[5], function(x) x/910),
lapply(Mb_table[6:7], function(x) x/880),
lapply(Mb_table[8:9], function(x) x/919),
lapply(Mb_table[10:16], function(x) x/880),
lapply(Mb_table[17:18], function(x) x/841), 
lapply(Mb_table[19:20], function(x) x/934),
lapply(Mb_table[21:22], function(x) x/851),
lapply(Mb_table[23], function(x) x/934),
lapply(Mb_table[28:29], function(x) x/1667),
lapply(Mb_table[30:31], function(x) x/1698),
lapply(Mb_table[32], function(x) x/1311) )

perc_table <- perc_table[,c(1:27)]
mydata <- t(perc_table)
scaledata <- scale(mydata, center = T)
d <- dist(scaledata, method = "euclidean")

cmdfit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(x, y, labels = row.names(mydata), cex=.7)
points(x, y, pch=19) 

cmdpoints <- as.data.frame(cmdfit$points)
cmdpoints$genome <- sub("_.*", "", row.names(cmdpoints))
cmdpoints$genome <- sub("ref", "", cmdpoints$genome)
cmdpoints$species <- rownames(cmdpoints)
cmdpoints$subsection <-cmdpoints$genome
cmdpoints$subsection <- gsub("D4|D7|D9|D11", "Erioxylum", cmdpoints$subsection)
cmdpoints$subsection <- gsub("D2.1|D2.2|D10", "Caducibracteata", cmdpoints$subsection)
cmdpoints$subsection <- gsub("D3D|D3K", "Integrifolia", cmdpoints$subsection)
cmdpoints$subsection <- gsub("D5", "Austroamericana", cmdpoints$subsection)
cmdpoints$subsection <- gsub("D6", "Selera", cmdpoints$subsection)
cmdpoints$subsection <- gsub("D1|D8", "Houzingenia", cmdpoints$subsection)
#cmdpoints$subsection <- gsub("A1|A2", "Gossypium", cmdpoints$subsection)
#cmdpoints$subsection <- gsub("F1", "Longiloba", cmdpoints$subsection)



# Erioxylum: D4, D7, D9, D11
# Caducibracteata: D2-1, D2-2, D10
# Integrifolia: D3d, D3k
# Austroamericana: D5
# Selera: D6
# Houzingenia: D1, D8
# Gossypium: A1, A2
# Longiloba: F1


# should add R hue here to differentiate colors
png("cotton.GS.ordination.png", 5000, 5000, pointsize=12, res=600)
ggplot(cmdpoints, aes(x=V1, y=V2, color=subsection)) + geom_point(size=1) + xlab("PCoA component 1") + ylab("PCoA component 2") + geom_text_repel(aes(label=species))+stat_ellipse()
dev.off()

### PCoA of log-transformed data ###
# evaluates how close the overall repetitive profiles are to one another

ord_table <- ord_table[,c(1:27)]
ord_table <- ord_table[!(rowSums(ord_table==0)),]

ord_table[ord_table==0]=0.00000001 # no -Inf values in PCA, so make them really small here
logdata <- t(log(ord_table))

dlog <- dist(logdata, method = "euclidean")

logcmd <- cmdscale(dlog,eig=TRUE, k=2) # k is the number of dim
logx <- logcmd$points[,1]
logy <- logcmd$points[,2]
plot(logx, logy, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(logx, logy, labels = row.names(logdata), cex=.7)
points(logx, logy, pch=19) 

logcmdpoints <- as.data.frame(logcmd$points)
logcmdpoints$genome <- sub("_.*", "", row.names(logcmdpoints))
logcmdpoints$genome <- sub("ref", "", logcmdpoints$genome)
logcmdpoints$species <- rownames(logcmdpoints)
logcmdpoints$subsection <-logcmdpoints$genome
logcmdpoints$subsection <- gsub("D4|D7|D9|D11", "Erioxylum", logcmdpoints$subsection)
logcmdpoints$subsection <- gsub("D2.1|D2.2|D10", "Caducibracteata", logcmdpoints$subsection)
logcmdpoints$subsection <- gsub("D3D|D3K", "Integrifolia", logcmdpoints$subsection)
logcmdpoints$subsection <- gsub("D5", "Austroamericana", logcmdpoints$subsection)
logcmdpoints$subsection <- gsub("D6", "Selera", logcmdpoints$subsection)
logcmdpoints$subsection <- gsub("D1|D8", "Houzingenia", logcmdpoints$subsection)

png("cotton.ordination.log.png", 5000, 5000, pointsize=12, res=600)
ggplot(logcmdpoints, aes(x=V1, y=V2, color=subsection)) + geom_point(size=2) + xlab("PCoA component 1") + ylab("PCoA component 2") + geom_text_repel(aes(label=species))+stat_ellipse()
dev.off()


### PCA to get the variation ###

cluster.pca <- prcomp(logdata, scale = TRUE)
cluster.eig <- get_eigenvalue(cluster.pca)

ind.coord <- cluster.pca$x

subsections <- as.data.frame(logdata[,1:2])
subsections$sub <-c("Houzingenia", "Houzingenia", 
"Caducibracteata", "Caducibracteata", 
"Integrifolia", "Integrifolia", "Integrifolia",
"Erioxylum", "Erioxylum",
"Austroamericana","Austroamericana","Austroamericana","Austroamericana","Austroamericana","Austroamericana","Austroamericana",
"Selera", "Selera", 
"Erioxylum","Erioxylum",
"Houzingenia", "Houzingenia", 
"Erioxylum",
"Caducibracteata","Caducibracteata","Caducibracteata",
"Erioxylum")

subfac <- as.factor(subsections[,3])

png("cotton_outgroup.PCA.direct.annot.png", 5000, 5000, pointsize=12, res=600)
fviz_pca_ind(cluster.pca, habillage=subfac, pointsize =2, invisible="quali", repel=TRUE, labelsize=3) + theme_minimal() + labs(title = "PCA of log counts") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text(face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5), legend.position="none") +theme_set(theme_grey(base_size=12))
dev.off()


### took out Procrustea ANOVA because nothing appears different
# maybe it should go back in to differentiate subsections?

### which species differ from one another in repeats
### https://groups.google.com/forum/#!topic/geomorph-r-package/8_B2thhxU_o

# Adams, D.C., and E. Otarola-Castillo. 2013. geomorph: an R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution. 4:393-399.
# Adams, D.C., M. L. Collyer, A. Kaliontzopoulou, and E. Sherratt. 2016 geomorph: Software for geometric morphometric analyses. R package version 3.0.2 http://cran.r-project.org/web/packages/geomorph/index.html.

# procD.lm with data as raw count data
mycountdata <- as.matrix(t(ord_table))
count.diff <- advanced.procD.lm(mycountdata~subfac, ~1, groups=~subfac)

#Effect sizes (Z)
#                Austroamericana Caducibracteata Erioxylum Houzingenia Integrifolia    Selera
#Austroamericana        0.000000       1.5480375 1.4458440   1.1730662    1.4247672 1.1961982
#Caducibracteata        1.548037       0.0000000 0.8888284   1.1820123    1.0312718 1.0868600
#Erioxylum              1.445844       0.8888284 0.0000000   1.2057821    0.9522423 1.1302263
#Houzingenia            1.173066       1.1820123 1.2057821   0.0000000    1.2384624 0.7678256
#Integrifolia           1.424767       1.0312718 0.9522423   1.2384624    0.0000000 1.2668033
#Selera                 1.196198       1.0868600 1.1302263   0.7678256    1.2668033 0.0000000

#P-values
#                Austroamericana Caducibracteata Erioxylum Houzingenia Integrifolia Selera
#Austroamericana           1.000           0.033     0.047       0.201        0.080  0.172
#Caducibracteata           0.033           1.000     0.565       0.211        0.333  0.247
#Erioxylum                 0.047           0.565     1.000       0.181        0.435  0.221
#Houzingenia               0.201           0.211     0.181       1.000        0.192  0.637
#Integrifolia              0.080           0.333     0.435       0.192        1.000  0.171
#Selera                    0.172           0.247     0.221       0.637        0.171  1.000


########### characterize composition ###########

annot_clust$cluster <- NULL

# 9.5 multiplier represents # reads (x) * 95nt/read * 1 kb/1000nt * 100% = # reads * 9.5 = # Kb in entire genome for that class 
Kbamount <- data.frame(annot_clust[1], apply(annot_clust[2:33], 2, function (x) x*9.5))
KBsum <- aggregate(. ~Lineage, data=Kbamount, FUN=sum)

KBsum$D1sum  <- rowMeans(KBsum[,2:3])
KBsum$D2.1sum <-KBsum$D2.1_6
KBsum$D2.2sum <-KBsum$D2.2
KBsum$D3Dsum <-KBsum$D3D_2
KBsum$D3Ksum <- rowMeans(KBsum[,7:8])
KBsum$D4sum  <- rowMeans(KBsum[,9:10])
KBsum$D5sum  <- rowMeans(KBsum[,11:17])
KBsum$D6sum  <- rowMeans(KBsum[,18:19])
KBsum$D7sum  <- rowMeans(KBsum[,20:21])
KBsum$D8sum  <- rowMeans(KBsum[,22:23])
KBsum$D9sum <-KBsum$D9_4
KBsum$D10sum <- rowMeans(KBsum[,25:27])
KBsum$D11sum <-KBsum$D11_1
KBsum$A1sum  <- rowMeans(KBsum[,29:30])
KBsum$A2sum  <- rowMeans(KBsum[,31:32])
KBsum$F1sum <-KBsum$F1_1

KBsum$D1min  <- apply(KBsum[,2:3], 1, min)
KBsum$D2.1min <-KBsum$D2.1_6
KBsum$D2.2min <-KBsum$D2.2
KBsum$D3Dmin <-KBsum$D3D_2
KBsum$D3Kmin <- apply(KBsum[,7:8], 1, min)
KBsum$D4min  <- apply(KBsum[,9:10], 1, min)
KBsum$D5min  <- apply(KBsum[,11:17], 1, min)
KBsum$D6min  <- apply(KBsum[,18:19], 1, min)
KBsum$D7min  <- apply(KBsum[,20:21], 1, min)
KBsum$D8min  <- apply(KBsum[,22:23], 1, min)
KBsum$D9min <-KBsum$D9_4
KBsum$D10min <- apply(KBsum[,25:27], 1, min)
KBsum$D11min <-KBsum$D11_1
KBsum$A1min  <- apply(KBsum[,29:30], 1, min)
KBsum$A2min  <- apply(KBsum[,31:32], 1, min)
KBsum$F1min <-KBsum$F1_1

KBsum$D1max  <- apply(KBsum[,2:3], 1, max)
KBsum$D2.1max <-KBsum$D2.1_6
KBsum$D2.2max <-KBsum$D2.2
KBsum$D3Dmax <-KBsum$D3D_2
KBsum$D3Kmax <- apply(KBsum[,7:8], 1, max)
KBsum$D4max  <- apply(KBsum[,9:10], 1, max)
KBsum$D5max  <- apply(KBsum[,11:17], 1, max)
KBsum$D6max  <- apply(KBsum[,18:19], 1, max)
KBsum$D7max  <- apply(KBsum[,20:21], 1, max)
KBsum$D8max  <- apply(KBsum[,22:23], 1, max)
KBsum$D9max <-KBsum$D9_4
KBsum$D10max <- apply(KBsum[,25:27], 1, max)
KBsum$D11max <-KBsum$D11_1
KBsum$A1max  <- apply(KBsum[,29:30], 1, max)
KBsum$A2max  <- apply(KBsum[,31:32], 1, max)
KBsum$F1max <-KBsum$F1_1

KBsum <- KBsum[,-(2:33)]
KBm <- melt(KBsum[,-(18:49)])
min <- c(KBsum$D1min,
KBsum$D2.1min,
KBsum$D2.2min,
KBsum$D3Dmin,
KBsum$D3Kmin,
KBsum$D4min,
KBsum$D5min,
KBsum$D6min,
KBsum$D7min,
KBsum$D8min,
KBsum$D9min,
KBsum$D10min,
KBsum$D11min,
KBsum$A1min,
KBsum$A2min,
KBsum$F1min)
max <- c(KBsum$D1max,
KBsum$D2.1max,
KBsum$D2.2max,
KBsum$D3Dmax,
KBsum$D3Kmax,
KBsum$D4max,
KBsum$D5max,
KBsum$D6max,
KBsum$D7max,
KBsum$D8max,
KBsum$D9max,
KBsum$D10max,
KBsum$D11max,
KBsum$A1max,
KBsum$A2max,
KBsum$F1max)
KBm$min <- min
KBm$max <- max

limits <- aes(ymax=KBm$max, ymin=KBm$min)

dodge <- position_dodge(width=0.9)

png("Figure_TE.amounts.png", 7500, 5000, pointsize=12, res=600)

ggplot(KBm, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + scale_y_log10(labels=comma) + geom_errorbar(limits, position = dodge) + labs(title = "Aggregate amounts in each species", x="Broad element category", y="Aggregate amount (mean) in kilobases") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text( face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5))+theme_set(theme_grey(base_size=12))

dev.off()

D1<-sum(KBsum$D1sum)/1000
D2_1<-sum(KBsum$D2.1sum)/1000
D2_2<-sum(KBsum$D2.2sum)/1000
D3d<-sum(KBsum$D3Dsum)/1000
D3k<-sum(KBsum$D3Ksum)/1000
D4<-sum(KBsum$D4sum)/1000
D5<-sum(KBsum$D5sum)/1000
D6<-sum(KBsum$D6sum)/1000
D7<-sum(KBsum$D7sum)/1000
D8<-sum(KBsum$D8sum)/1000
D9<-sum(KBsum$D9sum)/1000
D10<-sum(KBsum$D10sum)/1000
D11<-sum(KBsum$D11sum)/1000
A1<-sum(KBsum$A1sum)/1000
A2<-sum(KBsum$A2sum)/1000
F1<-sum(KBsum$F1sum)/1000

names <- c("D1","D2-1","D2-2","D3d", "D3k", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "A1", "A2", "F1")
TEamounts <- as.data.frame(names)
TEamounts$totalTE <-(c(D1,D2_1,D2_2, D3d, D3k, D4, D5, D6, D7, D8, D9, D10, D11, A1, A2, F1))
TEamounts$genomeSize <- (c(841,856,910,910,880,919,880,841,934,851,934,910,929,1667,1698,1311))
TEamounts$percent <- TEamounts$totalTE/TEamounts$genomeSize

gypsies <- KBsum[7,(1:17)]
row.names(gypsies) <- gypsies$Lineage
gypsies$Lineage <- NULL
TEamounts$gypsies <- t(gypsies/1000)
TEamounts$gypPerc <- TEamounts$gypsies/TEamounts$genomeSize


# min(TEamounts$percent)
# [1] 0.4451014
# max(TEamounts$percent[c(1:13)])
# [1] 0.5196145

### as % genome size

# from above
# Kbamount <- data.frame(annot_clust[1], apply(annot_clust[2:33], 2, function (x) x*9.5))

# convert to percent of genome, asks are the relative proportions the same
#genome	size	reads
#D01	841	98942
#D02.1	856	100706
#D02.2	910	107059
#D10	910	107059
#D11	929	109295
#D3D	910	107059
#D3K	880	103530
#D4	919	108118
#D5	880	103530
#D6	841	98942
#D7	934	109883
#D8	851	100118
#D9	934	109883
#A1	1667	196118
#A2	1698	199765
#F1	1311	154236
GSamount <- data.frame(Kbamount$Lineage, lapply(Kbamount[2:3], function(x) 100*x/841000), 
lapply(Kbamount[4], function(x) 100*x/856000), 
lapply(Kbamount[5], function(x) 100*x/910000),
lapply(Kbamount[6], function(x) 100*x/910000),
lapply(Kbamount[7:8], function(x) 100*x/880000),
lapply(Kbamount[9:10], function(x) 100*x/919000),
lapply(Kbamount[11:17], function(x) 100*x/880000),
lapply(Kbamount[18:19], function(x) 100*x/841000), 
lapply(Kbamount[20:21], function(x) 100*x/934000),
lapply(Kbamount[22:23], function(x) 100*x/851000),
lapply(Kbamount[24], function(x) 100*x/934000),
lapply(Kbamount[25:27], function(x) 100*x/910000),
lapply(Kbamount[28], function(x) 100*x/929000),lapply(Kbamount[29:30], function(x) 100*x/1667000),
lapply(Kbamount[31:32], function(x) 100*x/1698000),
lapply(Kbamount[33], function(x) 100*x/1311000) )

GSsum <- aggregate(. ~Kbamount.Lineage, data=GSamount, FUN=sum)

GSsum$D1sum  <- rowMeans(GSsum[,2:3])
GSsum$D2.1sum <-GSsum$D2.1_6
GSsum$D2.2sum <-GSsum$D2.2
GSsum$D3Dsum <-GSsum$D3D_2
GSsum$D3Ksum <- rowMeans(GSsum[,7:8])
GSsum$D4sum  <- rowMeans(GSsum[,9:10])
GSsum$D5sum  <- rowMeans(GSsum[,11:17])
GSsum$D6sum  <- rowMeans(GSsum[,18:19])
GSsum$D7sum  <- rowMeans(GSsum[,20:21])
GSsum$D8sum  <- rowMeans(GSsum[,22:23])
GSsum$D9sum <-GSsum$D9_4
GSsum$D10sum <- rowMeans(GSsum[,25:27])
GSsum$D11sum <-GSsum$D11_1
GSsum$A1sum  <- rowMeans(GSsum[,29:30])
GSsum$A2sum  <- rowMeans(GSsum[,31:32])
GSsum$F1sum <-GSsum$F1_1

GSsum$D1min  <- apply(GSsum[,2:3], 1, min)
GSsum$D2.1min <-GSsum$D2.1_6
GSsum$D2.2min <-GSsum$D2.2
GSsum$D3Dmin <-GSsum$D3D_2
GSsum$D3Kmin <- apply(GSsum[,7:8], 1, min)
GSsum$D4min  <- apply(GSsum[,9:10], 1, min)
GSsum$D5min  <- apply(GSsum[,11:17], 1, min)
GSsum$D6min  <- apply(GSsum[,18:19], 1, min)
GSsum$D7min  <- apply(GSsum[,20:21], 1, min)
GSsum$D8min  <- apply(GSsum[,22:23], 1, min)
GSsum$D9min <-GSsum$D9_4
GSsum$D10min <- apply(GSsum[,25:27], 1, min)
GSsum$D11min <-GSsum$D11_1
GSsum$A1min  <- apply(GSsum[,29:30], 1, min)
GSsum$A2min  <- apply(GSsum[,31:32], 1, min)
GSsum$F1min <-GSsum$F1_1

GSsum$D1max  <- apply(GSsum[,2:3], 1, max)
GSsum$D2.1max <-GSsum$D2.1_6
GSsum$D2.2max <-GSsum$D2.2
GSsum$D3Dmax <-GSsum$D3D_2
GSsum$D3Kmax <- apply(GSsum[,7:8], 1, max)
GSsum$D4max  <- apply(GSsum[,9:10], 1, max)
GSsum$D5max  <- apply(GSsum[,11:17], 1, max)
GSsum$D6max  <- apply(GSsum[,18:19], 1, max)
GSsum$D7max  <- apply(GSsum[,20:21], 1, max)
GSsum$D8max  <- apply(GSsum[,22:23], 1, max)
GSsum$D9max <-GSsum$D9_4
GSsum$D10max <- apply(GSsum[,25:27], 1, max)
GSsum$D11max <-GSsum$D11_1
GSsum$A1max  <- apply(GSsum[,29:30], 1, max)
GSsum$A2max  <- apply(GSsum[,31:32], 1, max)
GSsum$F1max <-GSsum$F1_1

GSsum <- GSsum[,-(2:33)]
GSm <- melt(GSsum[,-(18:49)])
min <- c(GSsum$D1min,
GSsum$D2.1min,
GSsum$D2.2min,
GSsum$D3Dmin,
GSsum$D3Kmin,
GSsum$D4min,
GSsum$D5min,
GSsum$D6min,
GSsum$D7min,
GSsum$D8min,
GSsum$D9min,
GSsum$D10min,
GSsum$D11min,
GSsum$A1min,
GSsum$A2min,
GSsum$F1min)
max <- c(GSsum$D1max,
GSsum$D2.1max,
GSsum$D2.2max,
GSsum$D3Dmax,
GSsum$D3Kmax,
GSsum$D4max,
GSsum$D5max,
GSsum$D6max,
GSsum$D7max,
GSsum$D8max,
GSsum$D9max,
GSsum$D10max,
GSsum$D11max,
GSsum$A1max,
GSsum$A2max,
GSsum$F1max)
GSm$min <- min
GSm$max <- max

limits <- aes(ymax=GSm$max, ymin=GSm$min)

dodge <- position_dodge(width=0.9)

png("Figure_TE.amounts.GS.png", 7500, 5000, pointsize=12, res=600)

ggplot(GSm, aes(x=Kbamount.Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + geom_errorbar(limits, position = dodge) + labs(title = "Percent Genome Size", x="Broad element category", y="Percent Genome Size") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text( face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5))+theme_set(theme_grey(base_size=12))

dev.off()



###################### glm/anova for cluster abundance
#
glm_clust <- annot_clust[,c(2:28)]
subsections <- subfac

data.in <- t(glm_clust)

P.scores <- rep(NA, ncol(data.in))
for (i in 1:ncol(data.in)){
  P.scores[i] <- anova(lm(data.in[,i]~subsections))$`Pr(>F)`[1]
}

glm_clust$p.value <- P.scores

glm_clust$p.valueBH <- p.adjust(glm_clust$p.value, method="BH")

# get rid of NA
glm_clust$p.valueBH[is.na(glm_clust$p.valueBH)] <- 99

row.names(glm_clust[(glm_clust$p.valueBH<0.05),])

[1] "CL0020" "CL0067" "CL0096" "CL0112" "CL0142" "CL0172" "CL0181" "CL0186" "CL0211" "CL0248" "CL0253" "CL0254" "CL0295" "CL0351" "CL0360"


############################
## D6 clusters significantly different from other Houzingenia

D6Hglm <- annot_clust[,c(2:28)]
D6Hsub <- t(D6Hglm[c(1:3),])
D6Hsub[,1]<-"NW"
D6Hsub[17,1]<-"Gg"
D6Hsub[18,1]<-"Gg"
D6Hsub <- D6Hsub[,1]

data.in <- t(D6Hglm)
D6H.P.scores <- rep(NA, ncol(data.in))
for (i in 1:ncol(data.in)){
  D6H.P.scores[i] <- anova(lm(data.in[,i]~D6Hsub))$`Pr(>F)`[1]
}

D6Hglm$p.value <- D6H.P.scores

D6Hglm$p.valueBH <- p.adjust(D6Hglm$p.value, method="BH")

# get rid of NA
D6Hglm$p.valueBH[is.na(D6Hglm$p.valueBH)] <- 99

row.names(D6Hglm[(D6Hglm$p.valueBH<0.05),])
# Notably, there are none

D6Hsig <- row.names(D6Hglm[(D6Hglm$p.valueBH<0.1),])
# "CL0186" "CL0368"

## D6 clusters significantly similar to Af species

D6Aglm <- annot_clust[,c(18:19,29:33)]
D6Asub <- t(D6Aglm[c(1:3),])
D6Asub[,1]<-"Af"
D6Asub[1,1]<-"Gg"
D6Asub[2,1]<-"Gg"
D6Asub <- D6Asub[,1]

data.in <- t(D6Aglm)
D6A.P.scores <- rep(NA, ncol(data.in))
for (i in 1:ncol(data.in)){
  D6A.P.scores[i] <- anova(lm(data.in[,i]~D6Asub))$`Pr(>F)`[1]
}

D6Aglm$p.value <- D6A.P.scores

D6Aglm$p.valueBH <- p.adjust(D6Aglm$p.value, method="BH")

# get rid of NA
D6Aglm$p.valueBH[is.na(D6Aglm$p.valueBH)] <- 99

row.names(D6Aglm[(D6Aglm$p.valueBH<0.05),])
#  [1] "CL0003" "CL0004" "CL0005" "CL0006" "CL0007" "CL0008" "CL0011" "CL0012" "CL0013" "CL0015" "CL0017" "CL0018" "CL0019" "CL0024" "CL0025" "CL0026" "CL0027" "CL0028" "CL0029" "CL0030" "CL0031" "CL0034" "CL0035" "CL0039"
# [25] "CL0041" "CL0042" "CL0043" "CL0046" "CL0049" "CL0052" "CL0054" "CL0055" "CL0056" "CL0059" "CL0060" "CL0062" "CL0064" "CL0067" "CL0069" "CL0071" "CL0073" "CL0074" "CL0075" "CL0077" "CL0079" "CL0082" "CL0083" "CL0084"
# [49] "CL0086" "CL0087" "CL0088" "CL0089" "CL0090" "CL0092" "CL0095" "CL0096" "CL0097" "CL0098" "CL0099" "CL0100" "CL0102" "CL0103" "CL0105" "CL0106" "CL0109" "CL0112" "CL0113" "CL0118" "CL0120" "CL0122" "CL0123" "CL0124"
# [73] "CL0127" "CL0128" "CL0129" "CL0134" "CL0136" "CL0140" "CL0142" "CL0145" "CL0148" "CL0149" "CL0151" "CL0152" "CL0154" "CL0155" "CL0156" "CL0158" "CL0159" "CL0160" "CL0161" "CL0162" "CL0163" "CL0164" "CL0171" "CL0173"
# [97] "CL0174" "CL0176" "CL0180" "CL0182" "CL0184" "CL0186" "CL0187" "CL0191" "CL0199" "CL0200" "CL0203" "CL0204" "CL0206" "CL0211" "CL0212" "CL0213" "CL0214" "CL0215" "CL0217" "CL0219" "CL0220" "CL0221" "CL0227" "CL0230"
#[121] "CL0231" "CL0248" "CL0250" "CL0261" "CL0277" "CL0278" "CL0280" "CL0284" "CL0285" "CL0286" "CL0287" "CL0288" "CL0293" "CL0296" "CL0297" "CL0301" "CL0302" "CL0303" "CL0305" "CL0306" "CL0308" "CL0313" "CL0314" "CL0319"
#[145] "CL0329" "CL0330" "CL0333" "CL0338" "CL0340" "CL0341" "CL0342" "CL0353" "CL0356" "CL0357" "CL0365" "CL0368" "CL0371" "CL0372" "CL0373" "CL0375" "CL0384" "CL0385" "CL0387" "CL0389" "CL0390" "CL0391"

# notably "CL0186" and "CL0368" also show up here

### Af versus Houz

D6AZglm <- annot_clust[,c(2:17,20:33)]
D6AZsub <- t(D6AZglm[c(1:3),])
D6AZsub[,1]<-"Hz"
D6AZsub[c(26:30),1]<-"Af"
D6AZsub <- D6AZsub[,1]

data.in <- t(D6AZglm)
D6AZ.P.scores <- rep(NA, ncol(data.in))
for (i in 1:ncol(data.in)){
  D6AZ.P.scores[i] <- anova(lm(data.in[,i]~D6AZsub))$`Pr(>F)`[1]
}

D6AZglm$p.value <- D6AZ.P.scores

D6AZglm$p.valueBH <- p.adjust(D6AZglm$p.value, method="BH")

# get rid of NA
D6AZglm$p.valueBH[is.na(D6AZglm$p.valueBH)] <- 99

row.names(D6AZglm[(D6AZglm$p.valueBH<0.05),])
#  [1] "CL0003" "CL0004" "CL0005" "CL0006" "CL0007" "CL0008" "CL0011" "CL0012" "CL0013" "CL0015" "CL0017" "CL0018" "CL0019" "CL0024" "CL0025" "CL0026" "CL0027" "CL0028" "CL0029" "CL0030" "CL0031" "CL0034" "CL0035" "CL0039"
# [25] "CL0041" "CL0042" "CL0043" "CL0046" "CL0049" "CL0052" "CL0054" "CL0055" "CL0056" "CL0059" "CL0060" "CL0062" "CL0064" "CL0067" "CL0069" "CL0071" "CL0073" "CL0074" "CL0075" "CL0077" "CL0079" "CL0082" "CL0083" "CL0084"
# [49] "CL0086" "CL0087" "CL0088" "CL0089" "CL0090" "CL0092" "CL0095" "CL0096" "CL0097" "CL0098" "CL0099" "CL0100" "CL0102" "CL0103" "CL0105" "CL0106" "CL0109" "CL0112" "CL0113" "CL0118" "CL0120" "CL0122" "CL0123" "CL0124"
# [73] "CL0127" "CL0128" "CL0129" "CL0134" "CL0136" "CL0140" "CL0142" "CL0145" "CL0148" "CL0149" "CL0151" "CL0152" "CL0154" "CL0155" "CL0156" "CL0158" "CL0159" "CL0160" "CL0161" "CL0162" "CL0163" "CL0164" "CL0171" "CL0173"
# [97] "CL0174" "CL0176" "CL0180" "CL0182" "CL0184" "CL0186" "CL0187" "CL0191" "CL0199" "CL0200" "CL0203" "CL0204" "CL0206" "CL0211" "CL0212" "CL0213" "CL0214" "CL0215" "CL0217" "CL0219" "CL0220" "CL0221" "CL0227" "CL0230"
#[121] "CL0231" "CL0248" "CL0250" "CL0261" "CL0277" "CL0278" "CL0280" "CL0284" "CL0285" "CL0286" "CL0287" "CL0288" "CL0293" "CL0296" "CL0297" "CL0301" "CL0302" "CL0303" "CL0305" "CL0306" "CL0308" "CL0313" "CL0314" "CL0319"
#[145] "CL0329" "CL0330" "CL0333" "CL0338" "CL0340" "CL0341" "CL0342" "CL0353" "CL0356" "CL0357" "CL0365" "CL0368" "CL0371" "CL0372" "CL0373" "CL0375" "CL0384" "CL0385" "CL0387" "CL0389" "CL0390" "CL0391"

# notably "CL0186" and "CL0368" also show up here


######## divergence times ######## 

dsTime <- read.csv("final_dNdS.csv", header=T, row.names=1)[,c(2,4,6)]
dsTime <- dsTime[ !grepl("D6-5", dsTime$qspecies),]
dsTime <- dsTime[ !grepl("D6-5", dsTime$sspecies),]
dsTime$sspecies <- gsub("D9-4|D4-185|D7-157|D11-1", "Erio", dsTime$sspecies)
dsTime$qspecies <- gsub("D9-4|D4-185|D7-157|D11-1", "Erio", dsTime$qspecies)

dsTime$sspecies <- gsub("D8-8|D2-2|D3-K-57|D3-D-27|D5-8|D2-1-6|D10-7|D1-35", "other", dsTime$sspecies)
dsTime$qspecies <- gsub("D8-8|D2-2|D3-K-57|D3-D-27|D5-8|D2-1-6|D10-7|D1-35", "other", dsTime$qspecies)

dsTime <- subset(dsTime, sspecies != qspecies)
dsTime <- subset(dsTime, dS<3)

quantile(dsTime$dS)
#  0%     25%    50%    75%   100% 
#0.0000 0.0138 0.0223 0.0371 2.8869 

# median time for Malvaceae, as per Grover GBE 2018 = 3.61E-09

erioNodeMin = 0.0138/(2*3.61e-09)/1000000 # 1.911357
erioNodeMax = 0.0371/(2*3.61e-09)/1000000 # 5.138504
erioNodeMedian = 0.0223/(2*3.61e-09)/1000000 #3.088643


trfn = "small.concat.raxml.nwk.tre"
tr = read.tree(trfn)

# find the node that corresponds to the Erioxylum/rest of Houzingenia divergence; node 17
plot(tr, cex = 0.5)
nodelabels()

# add dates to the tree
cal <- makeChronosCalib(tr, node=17, age.min=erioNodeMin, age.max=erioNodeMin)
chr <- chronos(tr, calibration = cal)
write.tree(chr, file="phylo.date.tre", digits=10)



######## ancestral state reconstruction of genome size for figure ######## 

cattree <- read.nexus("small.concat.raxml.tre")

# cattree$tip.label
# [1] "D6_5"   "D8_8"   "D1_35"  "D3K_57" "D3D_27" "D5_8"   "D2_1_6" "D10_7"  "D2_2"   "D11_1"  "D7_157" "D4_185" "D9_4"   "F1_1"

cattree$tip.label <- gsub("2_","2-",cattree$tip.label)
cattree$tip.label <- gsub("_.*","",cattree$tip.label)
cattree$tip.label <- gsub("K","k",cattree$tip.label)
cattree$tip.label <- gsub("3D","3d",cattree$tip.label)
#cattree$tip.label <- gsub("q","_",cattree$tip.label)

gsv <- TEamounts[c(1:13,16),3]
names(gsv) <- TEamounts[c(1:13,16),1]

name.check(cattree,gsv)

checkModel <- function (tree, df, outdf="name")
{
	BM <- fitContinuous(tree, df , model="BM")
	OU <- fitContinuous(tree, df , model="OU")
	LA <- fitContinuous(tree, df , model="lambda")
	KA <- fitContinuous(tree, df , model="kappa")
	DE <- fitContinuous(tree, df , model="delta")
	EB <- fitContinuous(tree, df , model="EB")
	WH <- fitContinuous(tree, df , model="white")
	TR <- fitContinuous(tree, df , model="trend")

	df <- as.matrix(data.frame(mods=c("BM", "OU", "LA", "KA", "DE", "EB", "WH", "TR"), lnL="NA", aic="NA"))
	mdf <- df
	mdf[1,2] <- BM$opt$lnL
	mdf[2,2] <- OU$opt$lnL
	mdf[3,2] <- LA$opt$lnL
	mdf[4,2] <- KA$opt$lnL
	mdf[5,2] <- DE$opt$lnL
	mdf[6,2] <- EB$opt$lnL
	mdf[7,2] <- WH$opt$lnL
	mdf[8,2] <- TR$opt$lnL
	mdf[1,3] <- BM$opt$aic
	mdf[2,3] <- OU$opt$aic
	mdf[3,3] <- LA$opt$aic
	mdf[4,3] <- KA$opt$aic
	mdf[5,3] <- DE$opt$aic
	mdf[6,3] <- EB$opt$aic
	mdf[7,3] <- WH$opt$aic
	mdf[8,3] <- TR$opt$aic

	aic_all <- c(BM$opt$aic,LA$opt$aic,KA$opt$aic,DE$opt$aic)
	names(aic_all) <- c("BM", "LA", "KA", "DE")
	aicres <- aicw(aic_all)

	mdf <- apply(mdf[,2:3],2,as.numeric)
	row.names(mdf)=df[,1]
	
	assign(outdf, mdf, envir=globalenv())

	mdf <- as.data.frame(mdf)
	mdiff <- max(mdf$lnL)-min(mdf$lnL)

	assign(paste0(outdf,"diff"), mdiff, envir=globalenv())
	assign(paste0(outdf,"aicw"), aicres, envir=globalenv())
}

checkModel(cattree, gsv, outdf="GSphylog")

# > GSphylog
#          lnL      aic
# BM -79.55543 163.1109
# OU -79.55544 165.1109
# LA -79.55543 165.1109
# KA -79.55543 165.1109
# DE -73.78319 153.5664
# EB -74.77727 155.5545
# WH -86.02424 176.0485
# TR -76.43052 158.8610

# WH <- fitContinuous(cattree, genomesizes , model="white")
# WH$opt$sigsq
# [1] 12724.69

rescale.cat <- rescale(cattree, model="white", WH$opt$siqsq)
# This tree looks unreasonably like a rake; discarding this model


makeGSstate <- function (tree, df, name="name")
{
    GSgradient <- contMap(tree, df, res=1000, plot=FALSE, lwd=0.5, fsize=1, sig=0)
    GSgradient$tree$tip.label <- gsub("F1",paste0("G.longicalx, ",round(df[["F1"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D11",paste0("G.schwendimanii, ",round(df[["D11"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D6",paste0("G.gossypioides, ",round(df[["D6"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10",paste0("G.turnerii, ",round(df[["D10"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1",paste0("G.thurberi, ",round(df[["D1"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8",paste0("G.trilobum, ",round(df[["D8"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3d",paste0("G.davidsonii, ",round(df[["D3d"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3k",paste0("G.klotzschianum, ",round(df[["D3k"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5",paste0("G.raimondii, ",round(df[["D5"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2-2",paste0("G.harknessii, ",round(df[["D2-2"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2-1",paste0("G.armourianum, ",round(df[["D2-1"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7",paste0("G.lobatum, ",round(df[["D7"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4",paste0("G.aridum, ",round(df[["D4"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D9",paste0("G.laxum, ",round(df[["D9"]]), " Mbp"),GSgradient$tree$tip.label)
    GSfit <- fastAnc(tree, df, vars=TRUE, CI=TRUE)
    gfit <- round(GSfit$ace)
    assign(paste0("G",name), GSgradient, envir=globalenv())
    assign(paste0("fit",name),gfit,envir=globalenv())
}

makeGSstateF1 <- function (tree, df, name="name")
{
    GSgradient <- contMap(tree, df, res=1000, plot=FALSE, lwd=0.5, fsize=1, sig=0)
    GSgradient$tree$tip.label <- gsub("D6",paste0("G.gossypioides, ",round(df[["D6"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D11",paste0("G.schwendimanii, ",round(df[["D11"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10",paste0("G.turnerii, ",round(df[["D10"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1",paste0("G.thurberi, ",round(df[["D1"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8",paste0("G.trilobum, ",round(df[["D8"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3d",paste0("G.davidsonii, ",round(df[["D3d"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3k",paste0("G.klotzschianum, ",round(df[["D3k"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5",paste0("G.raimondii, ",round(df[["D5"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2-2",paste0("G.harknessii, ",round(df[["D2-2"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2-1",paste0("G.armourianum, ",round(df[["D2-1"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7",paste0("G.lobatum, ",round(df[["D7"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4",paste0("G.aridum, ",round(df[["D4"]]), " Mbp"),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D9",paste0("G.laxum, ",round(df[["D9"]]), " Mbp"),GSgradient$tree$tip.label)
    GSfit <- fastAnc(tree, df, vars=TRUE, CI=TRUE)
    gfit <- round(GSfit$ace)
    assign(paste0("G",name), GSgradient, envir=globalenv())
    assign(paste0("fit",name),gfit,envir=globalenv())
}

makeGSstate(cattree, gsv, name="GS.gradient")

noF1 <- drop.tip(cattree, "F1")
plot(noF1)
gsvF <- gsv[-14]

makeGSstateF1(noF1, gsvF, name="GS.F1")

GGS.gradient$tree$tip.label <- gsub("Mbp.*", "Mbp                _____________________", GGS.gradient$tree$tip.label)
GGS.F1$tree$tip.label <- gsub("Mbp.*", "Mbp", GGS.F1$tree$tip.label)

time.labels=c(6.58,2.56,1.91,1.76,1.56,0.73,1.45,0.01,1.17,0.54,1.13,1.05,0.92)

png("Figure.ancGS.png", 15000, 7320, pointsize=12, res=600)

layout(matrix(c(1,2,0,0), 1, 2, byrow=FALSE))
#layout.show(2)

cw<-reorder(GGS.gradient$tree)
par(fg="transparent", bg="white")
plot(GGS.gradient, lwd=c(5,1), outline=FALSE, fsize=1) + nodelabels(time.labels, adj=c(-0.15,0.2), frame="none", cex=1.5)
par(fg="black", bg="white")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
text(rep(max(obj$xx[1:Ntip(cw)]),Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label), font=3,pos=4,cex=2,offset=0.3)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])), rep(obj$yy[i],2),lty="dotted")

cw<-reorder(GGS.F1$tree)
par(fg="transparent")
plot(GGS.F1, lwd=c(5,1), outline=FALSE, fsize=1, direction="leftwards") + nodelabels(round(fitGS.gradient[-1]), adj=c(1.2,0.2), frame="none", cex=1.5)
dev.off()

png("Figure.ancGSscale.png", 10000, 7320, pointsize=12, res=600)
plot(GGS.F1)
dev.off()

######## ancestral state reconstruction of "very" significant clusters ######## 

cptree <- read.nexus("concat.raxml.tre")

plot.phylo(cptree, cex =1, label.offset=0.0005, align.tip.label=TRUE)
nodelabels(cptree$node, adj=c(1.1,-0.2), frame="none", cex=1)

sigclusters <- row.names(glm_clust[(glm_clust$p.valueBH<0.05),])
ancTable <- Kbamount[sigclusters ,-c(1,29:32)]
names(ancTable) <- gsub("\\.", "_", names(ancTable))
names(ancTable)[names(ancTable) == 'D5ref'] <- 'D5_ref'
names(ancTable)[names(ancTable) == 'D3D_2'] <- 'D3D_27'
names(ancTable)[names(ancTable) == 'D3K_5'] <- 'D3K_57'
names(ancTable)[names(ancTable) == 'D4_12'] <- 'D4_12C'

anc20 <- t(ancTable["CL0020", ])
anc67 <- t(ancTable["CL0067", ])
anc96 <- t(ancTable["CL0096", ])
anc112 <- t(ancTable["CL0112", ])
anc142 <- t(ancTable["CL0142", ])
anc172 <- t(ancTable["CL0172", ])
anc181 <- t(ancTable["CL0181", ])
anc186 <- t(ancTable["CL0186", ])
anc211 <- t(ancTable["CL0211", ])
anc248 <- t(ancTable["CL0248", ])
anc253 <- t(ancTable["CL0253", ])
anc254 <- t(ancTable["CL0254", ])
anc295 <- t(ancTable["CL0295", ])
anc351 <- t(ancTable["CL0351", ])
anc360 <- t(ancTable["CL0360", ])

varNames <- grep("anc", ls(), value=TRUE)
varNames <- varNames[c(2:16)]

for (anc in varNames) {  
	obj <- get(anc)
	names(obj) <- row.names(obj)
	assign(anc, obj, envir=globalenv())
}


name.check(cptree,anc20) # check one table to make sure the names match
#[1] "OK"



checkModel(cptree, anc20, outdf="mdf20")
checkModel(cptree, anc67, outdf="mdf67")
checkModel(cptree, anc96, outdf="mdf96")
checkModel(cptree, anc112, outdf="mdf112")
checkModel(cptree, anc142, outdf="mdf142")
checkModel(cptree, anc172, outdf="mdf172")
checkModel(cptree, anc181, outdf="mdf181")
checkModel(cptree, anc186, outdf="mdf186")
checkModel(cptree, anc211, outdf="mdf211")
checkModel(cptree, anc248, outdf="mdf248")
checkModel(cptree, anc253, outdf="mdf253")
checkModel(cptree, anc254, outdf="mdf254")
checkModel(cptree, anc295, outdf="mdf295")
checkModel(cptree, anc351, outdf="mdf351")
checkModel(cptree, anc360, outdf="mdf360")

data.frame(mdf20[,1],mdf67[,1],mdf96[,1],mdf112[,1],mdf142[,1],mdf172[,1],mdf181[,1],mdf186[,1],mdf211[,1],mdf248[,1],mdf253[,1],mdf254[,1],mdf295[,1],mdf351[,1],mdf360[,1])


row.names(as.data.frame(which(mdf96[,2] == min(mdf96[,2]))))

# how does the best model compare to the BM module; i.e., using AIC, how much information is lost if we use BM instead of the preferred model
bestModel <- data.frame(cluster=rep(NA, 15), bestModel=rep("", 15), modelAIC=rep("", 15), BMAIC=rep("", 15), stringsAsFactors=FALSE)

bestModel [1,] <- list((deparse(substitute(mdf20))), (row.names(as.data.frame(which(mdf20[,2] == min(mdf20[,2]))))), min(mdf20[,2]), (mdf20[1,2]))
bestModel [2,] <- list((deparse(substitute(mdf67))), (row.names(as.data.frame(which(mdf67[,2] == min(mdf67[,2]))))), min(mdf67[,2]), (mdf67[1,2]))
bestModel [3,] <- list((deparse(substitute(mdf96))), (row.names(as.data.frame(which(mdf96[,2] == min(mdf96[,2]))))), min(mdf96[,2]), (mdf96[1,2]))
bestModel [4,] <- list((deparse(substitute(mdf112))), (row.names(as.data.frame(which(mdf112[,2] == min(mdf112[,2]))))), min(mdf112[,2]), (mdf112[1,2]))
bestModel [5,] <- list((deparse(substitute(mdf142))), (row.names(as.data.frame(which(mdf142[,2] == min(mdf142[,2]))))), min(mdf142[,2]), (mdf142[1,2]))
bestModel [6,] <- list((deparse(substitute(mdf172))), (row.names(as.data.frame(which(mdf172[,2] == min(mdf172[,2]))))), min(mdf172[,2]), (mdf172[1,2]))
bestModel [7,] <- list((deparse(substitute(mdf181))), (row.names(as.data.frame(which(mdf181[,2] == min(mdf181[,2]))))), min(mdf181[,2]), (mdf181[1,2]))
bestModel [8,] <- list((deparse(substitute(mdf186))), (row.names(as.data.frame(which(mdf186[,2] == min(mdf186[,2]))))), min(mdf186[,2]), (mdf186[1,2]))
bestModel [9,] <- list((deparse(substitute(mdf211))), (row.names(as.data.frame(which(mdf211[,2] == min(mdf211[,2]))))), min(mdf211[,2]), (mdf211[1,2]))
bestModel [10,] <- list((deparse(substitute(mdf248))), (row.names(as.data.frame(which(mdf248[,2] == min(mdf248[,2]))))), min(mdf248[,2]), (mdf248[1,2]))
bestModel [11,] <- list((deparse(substitute(mdf253))), (row.names(as.data.frame(which(mdf253[,2] == min(mdf253[,2]))))), min(mdf253[,2]), (mdf253[1,2]))
bestModel [12,] <- list((deparse(substitute(mdf254))), (row.names(as.data.frame(which(mdf254[,2] == min(mdf254[,2]))))), min(mdf254[,2]), (mdf254[1,2]))
bestModel [13,] <- list((deparse(substitute(mdf295))), (row.names(as.data.frame(which(mdf295[,2] == min(mdf295[,2]))))), min(mdf295[,2]), (mdf295[1,2]))
bestModel [14,] <- list((deparse(substitute(mdf351))), (row.names(as.data.frame(which(mdf351[,2] == min(mdf351[,2]))))), min(mdf351[,2]), (mdf351[1,2]))
bestModel [15,] <- list((deparse(substitute(mdf360))), (row.names(as.data.frame(which(mdf360[,2] == min(mdf360[,2]))))), min(mdf360[,2]), (mdf360[1,2]))

bestModel <- bestModel[order(bestModel$bestModel),]

bestModel$diff <- abs(as.numeric(bestModel$modelAIC) - as.numeric(bestModel$BMAIC))

# > bestModel 
#    cluster bestModel         modelAIC            BMAIC                        diff
# 15  mdf360        DE 387.460835793484 389.040430996053 0.4539366621825256520317282 # delta = 2.999999
# 3    mdf96        KA 444.731351462004  490.24071761278 0.0000000001311495192348412 # kappa = 0.000000
# 6   mdf172        KA 493.266114119696 515.140735620418 0.0000177822329084599169962 # kappa = 0.000000
# 7   mdf181        KA 394.868408744943 434.544847064236 0.0000000024230974924942434 # kappa = 0.000000
# 9   mdf211        KA  409.76815685806 419.786941026528 0.0066749598864533406933353 # kappa = 0.429722
# 13  mdf295        KA 440.691944049154 523.824343482667 0.0000000000000000008872149 # kappa = 0.000000
# 14  mdf351        KA 383.762863371115 454.250489449246 0.0000000000000004940904684 # kappa = 0.000000
# 1    mdf20        LA 406.751632340978 410.647522361433 0.1425667439648921064332399 # lambda = 0.994345 ~ Brownian motion
# 2    mdf67        LA 462.684227955132  512.98821550674 0.0000000000119296552758764 # lambda = 0.801979
# 4   mdf112        LA 436.906788912733 475.439564784961 0.0000000042925368600399812 # lambda = 0.896171
# 5   mdf142        LA 417.728198561328 457.495956007051 0.0000000023149477430454305 # lambda = 0.911297
# 8   mdf186        LA 428.100675567405 472.050171809862 0.0000000002860804298970051 # lambda = 0.826564
# 10  mdf248        LA 384.694961677806 398.501937705605 0.0010042763828229618258692 # lambda = 0.926967
# 11  mdf253        LA 393.279546824595 417.743195776489 0.0000048728845997401383453 # lambda = 0.949719
# 12  mdf254        LA 346.609100658292 410.092091343392 0.0000000000000163999953675 # lambda = 0.801065

# diff here represents the relative likelihood of each model
# https://en.wikipedia.org/wiki/Akaike_information_criterion
# As an example, suppose that there are three candidate models, whose AIC values are 100, 102, and 110. Then the second model is exp((100?-?102)/2) = 0.368 times 
# as probable as the first model to minimize the information loss. Similarly, the third model is exp((100?-?110)/2) = 0.007 times as probable as the first model to
# minimize the information loss.

# DE = delta is a time-dependent model of trait evolution (Pagel 1999). The delta model is similar
# to ACDC insofar as the delta model fits the relative contributions of early versus late evolution
# in the tree to the covariance of species trait values. Where delta is greater than 1, recent evolution
# has been relatively fast; if delta is less than 1, recent evolution has been comparatively
# slow. Intrepreted as a tree transformation, the model raises all node depths to an estimated
# power (delta). Default bounds are delta = c(min = exp(-500), max = 3)

# KA = kappa is a punctuational (speciational) model of trait evolution (Pagel 1999), where character
# divergence is related to the number of speciation events between two species. Note that if
# there are speciation events that are missing from the given phylogeny (due to extinction or
# incomplete sampling), interpretation under the kappa model may be difficult. Considered
# as a tree transformation, the model raises all branch lengths to an estimated power (kappa).
# Default bounds are kappa = c(min = exp(-500), max = 1) # {Pagel's kappa; raises all branch lengths to the power kappa. 
# As kappa approaches zero, the model becomes speciational.}

# LA = lambda is one of the Pagel (1999) models that fits the extent to which the phylogeny predicts
# covariance among trait values for species. The model effectively transforms the tree: values
# of lambda near 0 cause the phylogeny to become more star-like, and a lambda value of 1
# recovers the BM model. Default bounds are lambda = c(min = exp(-500), max = 1 # ? has a nice natural 
# scale between zero (no correlation between species) and 1.0 (correlation between species equal to the  
# Brownian expectation)

# > mdf20aicw
#         fit     delta          w
# BM 410.6475 3.8958900 0.06761687
# LA 406.7516 0.0000000 0.47428222
# KA 411.3453 4.5937063 0.04770086
# DE 407.0410 0.2893402 0.41040005
# > mdf67aicw
#         fit     delta                    w
# BM 512.9882 50.303988 0.000000000008432178
# LA 462.6842  0.000000 0.706824916469716125
# KA 464.4443  1.760026 0.293175083128668745
# DE 505.3038 42.619547 0.000000000393183033
# > mdf96aicw
#         fit       delta                   w
# BM 490.2407 45.50936615 0.00000000006605064
# LA 444.7604  0.02902899 0.49637143818390328
# KA 444.7314  0.00000000 0.50362855784002514
# DE 482.0790 37.34761366 0.00000000391002097
# > mdf112aicw
#         fit     delta                 w
# BM 475.4396 38.532776 0.000000003408803
# LA 436.9068  0.000000 0.794123091237706
# KA 439.6067  2.699921 0.205876822648245
# DE 469.0617 32.154932 0.000000082705247
# > mdf142aicw
#         fit     delta                 w
# BM 457.4960 39.767757 0.000000001774635
# LA 417.7282  0.000000 0.766598417127643
# KA 420.1066  2.378405 0.233401556771063
# DE 452.2600 34.531802 0.000000024326659
# > mdf172aicw
#         fit     delta             w
# BM 515.1407 21.874622 0.00001707236
# LA 499.6928  6.426713 0.03861571383
# KA 493.2661  0.000000 0.96007943898
# DE 506.4943 13.228200 0.00128777483
# > mdf181aicw
#         fit     delta                 w
# BM 434.5448 39.676438 0.000000002093477
# LA 398.5657  3.697279 0.136032686170675
# KA 394.8684  0.000000 0.863967201621429
# DE 426.6195 31.751051 0.000000110114420
# > mdf186aicw
#         fit     delta                 w
# BM 472.0502 43.949496 0.000000000232940
# LA 428.1007  0.000000 0.814246414539863
# KA 431.0564  2.955684 0.185753577405213
# DE 465.0223 36.921671 0.000000007821984
# > mdf211aicw
#         fit     delta           w
# BM 419.7869 10.018784 0.004766394
# LA 411.7525  1.984367 0.264753394
# KA 409.7682  0.000000 0.714070795
# DE 417.3144  7.546253 0.016409418
# > mdf248aicw
#         fit      delta            w
# BM 398.5019 13.8069760 0.0005629497
# LA 384.6950  0.0000000 0.5605526036
# KA 385.1863  0.4913263 0.4384562150
# DE 399.0490 14.3540285 0.0004282316
# > mdf253aicw
#         fit     delta              w
# BM 417.7432 24.463649 0.000003391733
# LA 393.2795  0.000000 0.696042133437
# KA 394.9372  1.657697 0.303858355228
# DE 411.0547 17.775144 0.000096119602
# > mdf254aicw
#         fit     delta                     w
# BM 410.0921 63.482991 0.0000000000000101236
# LA 346.6091  0.000000 0.6172926239846313612
# KA 347.5652  0.956145 0.3827073760143643866
# DE 400.9183 54.309179 0.0000000000009940375
# > mdf295aicw
#         fit     delta                          w
# BM 523.8243 83.132399 0.000000000000000000780959
# LA 444.6813  3.989345 0.119763426800766265012399
# KA 440.6919  0.000000 0.880236573199233540698572
# DE 513.9871 73.295117 0.000000000000000106848122
# > mdf351aicw
#         fit     delta                        w
# BM 454.2505 70.487626 0.0000000000000004905396
# LA 393.6195  9.856636 0.0071866452343565447730
# KA 383.7629  0.000000 0.9928133547655831891987
# DE 444.6448 60.881888 0.0000000000000597769226
# > mdf360aicw
#         fit    delta          w
# BM 389.0404 1.579595 0.25187263
# LA 391.0404 3.579595 0.09265876
# KA 390.8759 3.415025 0.10060567
# DE 387.4608 0.000000 0.55486294

###
# based on the above, BM seems reasonable for most if you are considering the exp-diff.  In most cases, save cluster 360 and 20, the favored model is not that much more likely to minimize information loss. 
# for cluster 20, the lambda factor is so close to 1 (which is BM anyway), that we may as use BM. The one I am most concerned about for that is cluster 360. Interestingly, by straight up diff, cluster 360 
# isn't terrible, but some of the rest are.  Let's see if we can get models into this stuff...


delta.anc360 <- fitContinuous(cptree, anc360 , model="delta")
delta.tree360 <- rescale(cptree, "delta", delta.anc360$opt$delta)

kappa.anc96 <- fitContinuous(cptree, anc96 , model="kappa")
kappa.tree96 <- rescale(cptree, "kappa", kappa.anc96$opt$kappa)

kappa.anc172 <- fitContinuous(cptree, anc172 , model="kappa")
kappa.tree172 <- rescale(cptree, "kappa", kappa.anc172$opt$kappa)

kappa.anc181 <- fitContinuous(cptree, anc181 , model="kappa")
kappa.tree181 <- rescale(cptree, "kappa", kappa.anc181$opt$kappa)

kappa.anc211 <- fitContinuous(cptree, anc211 , model="kappa")
kappa.tree211 <- rescale(cptree, "kappa", kappa.anc211$opt$kappa)

kappa.anc295 <- fitContinuous(cptree, anc295 , model="kappa")
kappa.tree295 <- rescale(cptree, "kappa", kappa.anc295$opt$kappa)

kappa.anc351 <- fitContinuous(cptree, anc351 , model="kappa")
kappa.tree351 <- rescale(cptree, "kappa", kappa.anc351$opt$kappa)

lambda.anc20 <- fitContinuous(cptree, anc20, model="lambda")
lambda.tree20 <- rescale(cptree, "lambda", lambda.anc20$opt$lambda)

lambda.anc67 <- fitContinuous(cptree, anc67, model="lambda")
lambda.tree67 <- rescale(cptree, "lambda", lambda.anc67$opt$lambda)

lambda.anc112 <- fitContinuous(cptree, anc112, model="lambda")
lambda.tree112 <- rescale(cptree, "lambda", lambda.anc112$opt$lambda)

lambda.anc142 <- fitContinuous(cptree, anc142, model="lambda")
lambda.tree142 <- rescale(cptree, "lambda", lambda.anc142$opt$lambda)

lambda.anc186 <- fitContinuous(cptree, anc186, model="lambda")
lambda.tree186 <- rescale(cptree, "lambda", lambda.anc186$opt$lambda)

lambda.anc248 <- fitContinuous(cptree, anc248, model="lambda")
lambda.tree248 <- rescale(cptree, "lambda", lambda.anc248$opt$lambda)

lambda.anc253 <- fitContinuous(cptree, anc253, model="lambda")
lambda.tree253 <- rescale(cptree, "lambda", lambda.anc253$opt$lambda)

lambda.anc254 <- fitContinuous(cptree, anc254, model="lambda")
lambda.tree254 <- rescale(cptree, "lambda", lambda.anc254$opt$lambda)


# 15  mdf360        DE 387.460835793484 389.040430996053 0.4539366621825256520317282 # delta = 2.999999

# 3    mdf96        KA 444.731351462004  490.24071761278 0.0000000001311495192348412 # kappa = 0.000000
# 6   mdf172        KA 493.266114119696 515.140735620418 0.0000177822329084599169962 # kappa = 0.000000
# 7   mdf181        KA 394.868408744943 434.544847064236 0.0000000024230974924942434 # kappa = 0.000000
# 9   mdf211        KA  409.76815685806 419.786941026528 0.0066749598864533406933353 # kappa = 0.429722
# 13  mdf295        KA 440.691944049154 523.824343482667 0.0000000000000000008872149 # kappa = 0.000000
# 14  mdf351        KA 383.762863371115 454.250489449246 0.0000000000000004940904684 # kappa = 0.000000

# 1    mdf20        LA 406.751632340978 410.647522361433 0.1425667439648921064332399 # lambda = 0.994345 ~ Brownian motion
# 2    mdf67        LA 462.684227955132  512.98821550674 0.0000000000119296552758764 # lambda = 0.801979
# 4   mdf112        LA 436.906788912733 475.439564784961 0.0000000042925368600399812 # lambda = 0.896171
# 5   mdf142        LA 417.728198561328 457.495956007051 0.0000000023149477430454305 # lambda = 0.911297
# 8   mdf186        LA 428.100675567405 472.050171809862 0.0000000002860804298970051 # lambda = 0.826564
# 10  mdf248        LA 384.694961677806 398.501937705605 0.0010042763828229618258692 # lambda = 0.926967
# 11  mdf253        LA 393.279546824595 417.743195776489 0.0000048728845997401383453 # lambda = 0.949719
# 12  mdf254        LA 346.609100658292 410.092091343392 0.0000000000000163999953675 # lambda = 0.801065


# now that we have transformed trees, we can run ancestral state reconstruction on these

makeState <- function (tree, df, name="name")
{
    GSgradient <- contMap(tree, df, res=1000, plot=FALSE, lwd=0.5, fsize=1, sig=0)
    GSgradient$tree$tip.label <- gsub("F1_1",paste0("G.longicalx, ",round(df[["F1_1"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D6_7",paste0("G.gossypioides 7, ",round(df[["D6_7"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D6_5",paste0("G.gossypioides 5, ",round(df[["D6_5"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1_35",paste0("G.thurberi 35, ",round(df[["D1_35"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1_2",paste0("G.thurberi 2, ",round(df[["D1_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8_9",paste0("G.trilobum 9, ",round(df[["D8_9"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8_8",paste0("G.trilobum 8, ",round(df[["D8_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3D_27",paste0("G.davidsonii 27, ",round(df[["D3D_27"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3K_57",paste0("G.klotzschianum 57, ",round(df[["D3K_57"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3K_56",paste0("G.klotzschianum 56, ",round(df[["D3K_56"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_6",paste0("G.raimondii 6, ",round(df[["D5_6"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_53",paste0("G.raimondii 53, ",round(df[["D5_53"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_31",paste0("G.raimondii 31, ",round(df[["D5_31"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_2",paste0("G.raimondii 2, ",round(df[["D5_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_4",paste0("G.raimondii 4, ",round(df[["D5_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_8",paste0("G.raimondii 8, ",round(df[["D5_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_ref",paste0("G.raimondii ref, ",round(df[["D5_ref"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_3",paste0("G.turnerii 3, ",round(df[["D10_3"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_7",paste0("G.turnerii 7, ",round(df[["D10_7"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_8",paste0("G.turnerii 8, ",round(df[["D10_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_2",paste0("G.harknessii, ",round(df[["D2_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_1_6",paste0("G.armourianum 6, ",round(df[["D2_1_6"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D11_1",paste0("G.schwendimanii 1, ",round(df[["D11_1"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_157",paste0("G.lobatum 157, ",round(df[["D7_157"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_4",paste0("G.lobatum 7, ",round(df[["D7_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4_185",paste0("G.aridum 185, ",round(df[["D4_185"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D9_4",paste0("G.laxum 4, ",round(df[["D9_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4_12C",paste0("G.aridum 12C, ",round(df[["D4_12C"]])),GSgradient$tree$tip.label)
    GSfit <- fastAnc(tree, df, vars=TRUE, CI=TRUE)
    gfit <- round(GSfit$ace)
    assign(paste0("G",name), GSgradient, envir=globalenv())
    assign(paste0("fit",name),gfit,envir=globalenv())
}

makeState(delta.tree360, anc360, "anc360")
makeState(kappa.tree96,  anc96, "anc96")
makeState(kappa.tree172, anc172, "anc172")
makeState(kappa.tree181, anc181, "anc181")
makeState(kappa.tree211, anc211, "anc211")
makeState(kappa.tree295, anc295, "anc295")
makeState(kappa.tree351, anc351, "anc351")
makeState(lambda.tree20, anc20, "anc20")
makeState(lambda.tree67, anc67, "anc67")
makeState(lambda.tree112, anc112, "anc112")
makeState(lambda.tree142, anc142, "anc142")
makeState(lambda.tree186, anc186, "anc186")
makeState(lambda.tree248, anc248, "anc248")
makeState(lambda.tree253, anc253, "anc253")
makeState(lambda.tree254, anc254, "anc254")


png("Figure_grid.anc.png", 10000, 7500, pointsize=12, res=600)
split.screen(figs=c(4,4))

screen(1)

par(fg="transparent", bg="white")
plot(Ganc360, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc360), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="white")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc360$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")

screen(2)
par(fg="transparent")
plot(Ganc96, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc96), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc96$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")

screen(3)
par(fg="transparent")
plot(Ganc172, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc172), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc172$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(4)
par(fg="transparent")
plot(Ganc181, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc181), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc181$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(5)
par(fg="transparent")
plot(Ganc211, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc211), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc211$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(6)
par(fg="transparent")
plot(Ganc295, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc295), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc295$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(7)
par(fg="transparent")
plot(Ganc351, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc351), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc351$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(8)
par(fg="transparent")
plot(Ganc20, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc20), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc20$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(9)
par(fg="transparent")
plot(Ganc67, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc67), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc67$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(10)
par(fg="transparent")
plot(Ganc112, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc112), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc112$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(11)
par(fg="transparent")
plot(Ganc142, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc142), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc142$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(12)
par(fg="transparent")
plot(Ganc186, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc186), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc186$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(13)
par(fg="transparent")
plot(Ganc248, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc248), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc248$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(14)
par(fg="transparent")
plot(Ganc253, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc253), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc253$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(15)
par(fg="transparent")
plot(Ganc254, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitanc254), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cw<-reorder(Ganc254$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")

close.screen(all.screens=TRUE)

dev.off()








               
######## ancestral state reconstruction of copia clusters ######## 

cptree <- read.nexus("concat.raxml.tre")
cptree <- drop.tip(cptree, "'D10_8.o'")

#plot.phylo(cptree, cex =1, label.offset=0.0005, align.tip.label=TRUE)
#nodelabels(cptree$node, adj=c(1.1,-0.2), frame="none", cex=1)

copiaclusters <- row.names(annot_clust[(annot_clust$Lineage=="LTR/Copia"),])
copTable <- Kbamount[copiaclusters ,-c(1,29:32)]
names(copTable)[names(copTable) == 'D5ref'] <- 'D5_ref'
names(copTable)[names(copTable) == 'D3D_2'] <- 'D3D_27'
names(copTable)[names(copTable) == 'D3K_5'] <- 'D3K_57'
names(copTable)[names(copTable) == 'D4_12'] <- 'D4_12C'
names(copTable)[names(copTable) == 'D2.1_6'] <- 'D2_1_6'
names(copTable)[names(copTable) == 'D2.2'] <- 'D2_2'

#> copiaclusters
# [1] "CL0059" "CL0068" "CL0072" "CL0088" "CL0094" "CL0103" "CL0105" "CL0107" "CL0115" "CL0118"
#[11] "CL0125" "CL0126" "CL0129" "CL0132" "CL0139" "CL0141" "CL0144" "CL0146" "CL0150" "CL0159"
#[21] "CL0171" "CL0172" "CL0184" "CL0192" "CL0194" "CL0195" "CL0202" "CL0205" "CL0207" "CL0208"
#[31] "CL0219" "CL0227" "CL0250" "CL0253" "CL0257" "CL0271" "CL0283" "CL0294" "CL0295" "CL0298"
#[41] "CL0305" "CL0314" "CL0319" "CL0322" "CL0335" "CL0342" "CL0348" "CL0356" "CL0358" "CL0364"
#[51] "CL0382" "CL0388"

for (i in c(1:52)) { assign(paste0("cop",copiaclusters[i]), t(subset(copTable, row.names(copTable) %in% copiaclusters[i]))) }


varNames <- grep("copCL", ls(), value=TRUE)
#varNames <- varNames[c(2:16)]

for (cop in varNames) {  
	obj <- get(cop)
	names(obj) <- row.names(obj)
	assign(cop, obj, envir=globalenv())
}


name.check(cptree,copCL0059) # check one table to make sure the names match
#[1] "OK"


for (name in varNames) { checkModel(cptree, get(name),outdf=paste0("mdf",name)) }

mdfNames <- grep("mdfcopCL", ls(), value=TRUE)
mdfNames <- gsub("aicw|diff", "", mdfNames)
mdfNames <- mdfNames[!duplicated(mdfNames)]


# how does the best model compare to the BM module; i.e., using AIC, how much information is lost if we use BM instead of the preferred model
bestModel <- data.frame(cluster=rep(NA, length(mdfNames)), bestModel=rep("", length(mdfNames)), modelAIC=rep("", length(mdfNames)), BMAIC=rep("", length(mdfNames)), stringsAsFactors=FALSE)

for (i in (1:length(mdfNames))) {
bestModel[i,] <- list((mdfNames[i]), (row.names(as.data.frame(which(get(mdfNames[i])[,2] == min(get(mdfNames[i])[,2]))))), min(get(mdfNames[i])[,2]), (get(mdfNames[i])[1,2]))
 }

bestModel <- bestModel[order(bestModel$bestModel),]

bestModel$diff <- abs(as.numeric(bestModel$modelAIC) - as.numeric(bestModel$BMAIC))
bestModel$exp <- exp(bestModel$diff/2)

# Interesting that the WH model is most often the best, but let's not use it

makeState <- function (tree, df, name="name")
{
    GSgradient <- contMap(tree, df, res=1000, plot=FALSE, lwd=0.5, fsize=1, sig=0)
    GSgradient$tree$tip.label <- gsub("F1_1",paste0("G.longicalx, ",round(df[["F1_1"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D6_7",paste0("G.gossypioides 7, ",round(df[["D6_7"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D6_5",paste0("G.gossypioides 5, ",round(df[["D6_5"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1_35",paste0("G.thurberi 35, ",round(df[["D1_35"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1_2",paste0("G.thurberi 2, ",round(df[["D1_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8_9",paste0("G.trilobum 9, ",round(df[["D8_9"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8_8",paste0("G.trilobum 8, ",round(df[["D8_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3D_27",paste0("G.davidsonii 27, ",round(df[["D3D_27"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3K_57",paste0("G.klotzschianum 57, ",round(df[["D3K_57"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3K_56",paste0("G.klotzschianum 56, ",round(df[["D3K_56"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_6",paste0("G.raimondii 6, ",round(df[["D5_6"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_53",paste0("G.raimondii 53, ",round(df[["D5_53"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_31",paste0("G.raimondii 31, ",round(df[["D5_31"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_2",paste0("G.raimondii 2, ",round(df[["D5_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_4",paste0("G.raimondii 4, ",round(df[["D5_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_8",paste0("G.raimondii 8, ",round(df[["D5_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_ref",paste0("G.raimondii ref, ",round(df[["D5_ref"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_3",paste0("G.turnerii 3, ",round(df[["D10_3"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_7",paste0("G.turnerii 7, ",round(df[["D10_7"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_8",paste0("G.turnerii 8, ",round(df[["D10_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_2",paste0("G.harknessii, ",round(df[["D2_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_1_6",paste0("G.armourianum 6, ",round(df[["D2_1_6"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D11_1",paste0("G.schwendimanii 1, ",round(df[["D11_1"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_157",paste0("G.lobatum 157, ",round(df[["D7_157"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_4",paste0("G.lobatum 7, ",round(df[["D7_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4_185",paste0("G.aridum 185, ",round(df[["D4_185"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D9_4",paste0("G.laxum 4, ",round(df[["D9_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4_12C",paste0("G.aridum 12C, ",round(df[["D4_12C"]])),GSgradient$tree$tip.label)
    GSfit <- fastAnc(tree, df, vars=TRUE, CI=TRUE)
    gfit <- round(GSfit$ace)
    assign(paste0("G",name), GSgradient, envir=globalenv())
    assign(paste0("fit",name),gfit,envir=globalenv())
}

varNames <- gsub("mdf|aicw|diff", "", varNames)
varNames <- varNames [!duplicated(varNames )]

for (name in varNames) { makeState(cptree, get(name),(paste0("anc",name))) }

nodeSize <- ls(pattern="fitanc")

vals <- vector(mode="numeric", length=0) # compare tips to reconstructed basal D node
Fvals <- vector(mode="numeric", length=0)  # compare tips to reconstructed basal node
Fchange <- vector(mode="numeric", length=0)

for (name in row.names(copCL0059)) { 
    newvec <- vector(mode="numeric", length=0)
    assign(paste0("vec",name),newvec,envir=globalenv()) }

for (i in varNames) { 
    
    temp <- get(i)
    sub <- get(paste0("fitanc",i))

    # how much did each change compared to the reconstructed D root
    temp <- as.data.frame(temp[c(1:27),])
    temp$diff <- temp[,1] - sub[[2]]
    vals <-c(vals,temp$diff)    

    # how much did each change compared to reconstructed root
    temp$Fdiff <- temp[,1] - sub[[1]]
    Fvals <-c(Fvals,temp$Fdiff)

    for (r in c(1:27)) {
        newr <- temp[r,1]-sub[[2]]
        assign(paste0("vec",row.names(temp)[r]), c(get(paste0("vec",row.names(temp)[r])),newr), envir=globalenv()) 
    }   

    assign(i,temp,envir=globalenv())
}

vecs <- ls(pattern="vecD")
#  [1] "vecD1_2"   "vecD1_35"  "vecD10_3"  "vecD10_7"  "vecD10_8"  "vecD11_1"  "vecD2_1_6" "vecD2_2"   "vecD3D_27" "vecD3K_56" "vecD3K_57"
# [12] "vecD4_12C" "vecD4_185" "vecD5_2"   "vecD5_31"  "vecD5_4"   "vecD5_53"  "vecD5_6"   "vecD5_8"   "vecD5_ref" "vecD6_5"   "vecD6_7"  
# [23] "vecD7_157" "vecD7_4"   "vecD8_8"   "vecD8_9"   "vecD9_4"


boxplot(vals)
boxplot(Fvals)
boxplot(Fchange)

boxplot.matrix <- cbind(vals, Fvals,vecD1_2, vecD1_35, vecD2_1_6, vecD2_2, vecD3D_27, vecD3K_56, vecD3K_57, vecD4_12C, vecD4_185, vecD5_2, vecD5_31, vecD5_4, vecD5_53, vecD5_6, vecD5_8, vecD5_ref, vecD6_5, vecD6_7, vecD7_157, vecD7_4, vecD8_8, vecD8_9, vecD9_4, vecD10_3, vecD10_7, vecD10_8, vecD11_1, Fchange)


plusminus <- data.frame(underZero=numeric(), total_under=numeric(), overZero=numeric(), total_over=numeric(), median=numeric(), sum=numeric(), stringsAsFactors=FALSE)

plusminus[1,1] <- length(which(vals<0)) 
plusminus[1,2] <- sum(which(vals<0))
plusminus[1,3] <- length(which(vals>0)) 
plusminus[1,4] <- sum(which(vals>0))
plusminus[1,5] <- median(vals)
row.names(plusminus)[1] <- "overall"

for (i in vecs) {
    temp <- get(i)
    plusminus[i,1] <- length(temp[which(temp<0)]) 
    plusminus[i,2] <- sum(temp[which(temp<0)])
    plusminus[i,3] <- length(temp[which(temp>0)]) 
    plusminus[i,4] <- sum(temp[which(temp>0)])
    plusminus[i,5] <- median(temp)
    plusminus[i,6] <- sum(temp)
}

write.table(plusminus, file="copia.gainloss.tbl", quote=F, sep="\t", row.names=TRUE)

length(which(Fvals<0))
#[1] 664

length(which(Fvals>0))
#[1] 740

median(Fvals)
#[1] 12


gainloss <- c(sum(vecD1_2), sum(vecD1_35), sum(vecD2_1_6), sum(vecD2_2), sum(vecD3D_27), sum(vecD3K_56), sum(vecD3K_57), sum(vecD4_12C), sum(vecD4_185), sum(vecD5_2), sum(vecD5_31), sum(vecD5_4), sum(vecD5_53), sum(vecD5_6), sum(vecD5_8), sum(vecD5_ref), sum(vecD6_5), sum(vecD6_7), sum(vecD7_157), sum(vecD7_4), sum(vecD8_8), sum(vecD8_9), sum(vecD9_4), sum(vecD10_3), sum(vecD10_7), sum(vecD10_8), sum(vecD11_1))
names(gainloss)<- c("D1_2", "D1_35", "D2_1_6", "D2_2", "D3D_27", "D3K_56", "D3K_57", "D4_12C", "D4_185", "D5_2", "D5_31", "D5_4", "D5_53", "D5_6", "D5_8", "D5_ref", "D6_5", "D6_7", "D7_157", "D7_4", "D8_8", "D8_9", "D9_4", "D10_3", "D10_7", "D10_8", "D11_1")

#   D1_2   D1_35  D2_1_6    D2_2  D3D_27  D3K_56  D3K_57  D4_12C  D4_185    D5_2 
#-3995.5 -4765.0   650.0 -1962.5  1448.0  5694.5  8117.0  1913.5   298.5  3946.5 
#  D5_31    D5_4   D5_53    D5_6    D5_8  D5_ref    D6_5    D6_7  D7_157    D7_4 
#-2124.0 -3672.5 -5335.0 40379.0 -4841.0 -4470.5 -1145.5 -5012.0   545.5  1391.0 
#   D8_8    D8_9    D9_4   D10_3   D10_7   D10_8   D11_1 
# 3832.5 -3292.5  4070.0  -936.5 -4717.5 -3473.0 -1164.5 

gainloss <- melt(gainloss)

png("Figure.magnitude.gainloss.copia.png", 15000, 7320, pointsize=12, res=600)
ggplot(gainloss, aes(x=row.names(gainloss),y=value)) + geom_bar(stat="identity")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

png("Figure.number.gainloss.png", 15000, 7320, pointsize=12, res=600)
gldf <- as.data.frame(cbind(vals, Fvals, Fchange))
ggplot(gldf, aes(x=vals)) + geom_density(fill="blue",color="blue")+geom_vline(xintercept=0, size=1, linetype="dashed")
deorv.off()                             
                             


########### relative aging of transposable elements ###########

# see TE_dating_histogram.pl for table generation

ages <- read.table("TE_dating.txt", header=T, row.names=1, sep="\t")
ages$age <- ages$Category

allages <- read.table("dating.forR", header =T, row.names=1, sep="\t")
allages$age <- allages$Category

#Category	age
#1	young
#3	old
#4	young
#5	old
#4*	old

ages$age <- sub("3|5|4\\*","old",ages$age)
ages$age <- sub("1|4","young",ages$age)

allages$age <- sub("3|5|4\\*","old",allages$age)
allages$age <- sub("1|4","young",allages$age)

Dallages <- allages[(allages$name != "A1-155"),]
Dallages <- Dallages[(Dallages$name != "A1-97"),]
Dallages <- Dallages[(Dallages$name != "A2_101"),]
Dallages <- Dallages[(Dallages$name != "A2_34"),]
Dallages <- Dallages[(Dallages$name != "all"),]
Dallages <- Dallages[(Dallages$name != "F1-1"),]
Dallages <- Dallages[!is.na(Dallages$age),]


table(ages$age)
#  old young 
#  308    86 

> table(Dallages$age)

    2     6   old young 
    1    15  7979  2221 


table(subset(Dallages, name=="D4-185")$age)
#  old young 
#  286    79 


yold <- data.frame(old=numeric(), young=numeric(), perc_old=numeric(), stringsAsFactors=FALSE)

for (aname in Dnames) {
#    temp <- get(i)
    yold[aname,1] <- nrow(subset(Dallages, name==aname & age=="old")) 
    yold[aname,2] <- nrow(subset(Dallages, name==aname & age=="young"))
    yold[aname,3] <- 100*yold[aname,1]/(yold[aname,2] + yold[aname,1])
}

min(yold$perc_old)
# 76.5%
max(yold$perc_old)
# 79.7%

write.table(yold, file="D.ages.tbl", quote=F, sep="\t", row.names=TRUE)

DsmAge <- Dallages[,c(1,2,8)]
DsmAge <- recast(DsmAge, cluster~name)
DsmYoung <- DsmAge[rowSums(DsmAge == "young", na.rm = TRUE)>0,]
DsmYoung$number <- rowSums(DsmYoung == "young", na.rm = TRUE)
table(DsmYoung$number)
#  1  2  3  4  5  6  7  9 10 12 13 14 16 17 18 22 24 25 26 28 
# 20  5  7  2  3  1  4  1  3  1  1  3  1  1  1  1  2  4  3 61 

nrow(DsmYoung)
# 125

61/125
# 0.488
20/125
# 0.16

png("Figure.hist.number.young.png", 15000, 7320, pointsize=12, res=600)
hist(DsmYoung$number)
dev.off()

youngClusters <- as.character(DsmYoung$cluster)

youngTable <- Kbamount[youngClusters ,-c(1,29:32)]
names(youngTable)[names(youngTable) == 'D5ref'] <- 'D5_ref'
names(youngTable)[names(youngTable) == 'D3D_2'] <- 'D3D_27'
names(youngTable)[names(youngTable) == 'D3K_5'] <- 'D3K_57'
names(youngTable)[names(youngTable) == 'D4_12'] <- 'D4_12C'
names(youngTable)[names(youngTable) == 'D2.1_6'] <- 'D2_1_6'
names(youngTable)[names(youngTable) == 'D2.2'] <- 'D2_2'

for (i in c(1:125)) { assign(paste0("yng",youngClusters[i]), t(subset(youngTable, row.names(youngTable) %in% youngClusters[i]))) }


youngNames <- grep("yngCL", ls(), value=TRUE)


for (yng in youngNames) {  
	obj <- get(yng)
	names(obj) <- row.names(obj)
	assign(yng, obj, envir=globalenv())
}


name.check(cptree,yngCL0055) # check one table to make sure the names match
#[1] "OK"

for (name in youngNames) { makeState(cptree, get(name),(paste0("anc",name))) }

youngVals <- vector(mode="numeric", length=0) # compare tips to reconstructed basal D node

for (name in row.names(yngCL0055)) { 
    youngvec <- vector(mode="numeric", length=0)
    assign(paste0("young",name),newvec,envir=globalenv()) }

for (i in youngNames) { 
    
    temp <- get(i)
    sub <- get(paste0("fitanc",i))

    # how much did each change compared to the reconstructed D-F root
    temp <- as.data.frame(temp[c(1:27),])
    temp$diff <- temp[,1] - sub[[1]]
    youngVals <-c(youngVals,temp$diff)    

    for (r in c(1:27)) {
        newr <- temp[r,1]-sub[[1]]
        assign(paste0("young",row.names(temp)[r]), c(get(paste0("young",row.names(temp)[r])),newr), envir=globalenv()) 
    }   

    assign(i,temp,envir=globalenv())
}

youngs <- ls(pattern="youngD")

boxplot(youngVals)

yplusminus <- data.frame(underZero=numeric(), total_under=numeric(), overZero=numeric(), total_over=numeric(), median=numeric(), sum=numeric(), stringsAsFactors=FALSE)

for (i in youngs) {
    temp <- get(i)
    yplusminus[i,1] <- length(temp[which(temp<0)]) 
    yplusminus[i,2] <- sum(temp[which(temp<0)])
    yplusminus[i,3] <- length(temp[which(temp>0)]) 
    yplusminus[i,4] <- sum(temp[which(temp>0)])
    yplusminus[i,5] <- median(temp)
    yplusminus[i,6] <- sum(temp)
}

write.table(yplusminus, file="young.gainloss.total.tbl", quote=F, sep="\t", row.names=TRUE)

youngGrow <- DsmYoung[DsmYoung$number==28,]
growNames <- gsub("CL", "yngCL", youngGrow$cluster)


growVals <- vector(mode="numeric", length=0) # compare tips to reconstructed basal D-F node

for (i in growNames) { 
    
    temp <- get(i)
    sub <- get(paste0("fitanc",i))

    # how much did each change compared to the reconstructed D-F root
    temp <- as.data.frame(temp[c(1:27),])
    temp$diff <- temp[,1] - sub[[1]]
    growVals <-c(growVals,temp$diff)    

    for (r in c(1:27)) {
        newr <- temp[r,1]-sub[[1]]
        assign(paste0("grow",row.names(temp)[r]), c(get(paste0("young",row.names(temp)[r])),newr), envir=globalenv()) 
    }   

    assign(i,temp,envir=globalenv())
}

grows <- ls(pattern="growD")

boxplot(youngVals, growVals)

growplusminus <- data.frame(underZero=numeric(), total_under=numeric(), overZero=numeric(), total_over=numeric(), median=numeric(), sum=numeric(), stringsAsFactors=FALSE)

for (i in grows) {
    temp <- get(i)
    growplusminus[i,1] <- length(temp[which(temp<0)]) 
    growplusminus[i,2] <- sum(temp[which(temp<0)])
    growplusminus[i,3] <- length(temp[which(temp>0)]) 
    growplusminus[i,4] <- sum(temp[which(temp>0)])
    growplusminus[i,5] <- median(temp)
    growplusminus[i,6] <- sum(temp)
}


OneGrow <- DsmYoung[DsmYoung$number==1,]
OneNames <- gsub("CL", "yngCL", OneGrow$cluster)


Dclust <- annot_clust[,c(2:28)]
Dclust$sum <- rowSums(Dclust)
Dclust <- Dclust[Dclust$sum>279,]




makeStateFig <- function (tree, df, name="name")
{
    GSgradient <- contMap(tree, df, res=1000, plot=FALSE, lwd=0.5, fsize=1, sig=0)
    GSgradient$tree$tip.label <- gsub("F1_1",paste0("G.longicalx, ",round(df[["F1_1"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D6_7",paste0("G.gossypioides 7, ",round(df[["D6_7"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D6_5",paste0("G.gossypioides 5, ",round(df[["D6_5"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1_35",paste0("G.thurberi 35, ",round(df[["D1_35"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D1_2",paste0("G.thurberi 2, ",round(df[["D1_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8_9",paste0("G.trilobum 9, ",round(df[["D8_9"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D8_8",paste0("G.trilobum 8, ",round(df[["D8_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3D_27",paste0("G.davidsonii 27, ",round(df[["D3D_27"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3K_57",paste0("G.klotzschianum 57, ",round(df[["D3K_57"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D3K_56",paste0("G.klotzschianum 56, ",round(df[["D3K_56"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_6",paste0("G.raimondii 6, ",round(df[["D5_6"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_53",paste0("G.raimondii 53, ",round(df[["D5_53"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_31",paste0("G.raimondii 31, ",round(df[["D5_31"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_2",paste0("G.raimondii 2, ",round(df[["D5_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_4",paste0("G.raimondii 4, ",round(df[["D5_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_8",paste0("G.raimondii 8, ",round(df[["D5_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D5_ref",paste0("G.raimondii ref, ",round(df[["D5_ref"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_3",paste0("G.turnerii 3, ",round(df[["D10_3"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_7",paste0("G.turnerii 7, ",round(df[["D10_7"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_8",paste0("G.turnerii 8, ",round(df[["D10_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_2",paste0("G.harknessii, ",round(df[["D2_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_1_6",paste0("G.armourianum 6, ",round(df[["D2_1_6"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D11_1",paste0("G.schwendimanii 1, ",round(df[["D11_1"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_157",paste0("G.lobatum 157, ",round(df[["D7_157"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_4",paste0("G.lobatum 7, ",round(df[["D7_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4_185",paste0("G.aridum 185, ",round(df[["D4_185"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D9_4",paste0("G.laxum 4, ",round(df[["D9_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4_12C",paste0("G.aridum 12C, ",round(df[["D4_12C"]])),GSgradient$tree$tip.label)
    GSfit <- fastAnc(tree, df, vars=TRUE, CI=TRUE)
    gfit <- round(GSfit$ace)
    
    png(paste0(name,".anc.png"), 5000, 3750, pointsize=12, res=600)
    plot(GSgradient, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(gfit), adj=c(-0.3,0.2), frame="none", cex=0.6)
    dev.off()

    assign(paste0("G",name), GSgradient, envir=globalenv())
    assign(paste0("fit",name),gfit,envir=globalenv())
}

allTable <- Kbamount[,-c(1,29:32)]
names(allTable)[names(allTable) == 'D5ref'] <- 'D5_ref'
names(allTable)[names(allTable) == 'D3D_2'] <- 'D3D_27'
names(allTable)[names(allTable) == 'D3K_5'] <- 'D3K_57'
names(allTable)[names(allTable) == 'D4_12'] <- 'D4_12C'
names(allTable)[names(allTable) == 'D2.1_6'] <- 'D2_1_6'
names(allTable)[names(allTable) == 'D2.2'] <- 'D2_2'


for (i in c(1:nrow(allTable))) { assign(paste0("all",row.names(allTable)[i]), t(allTable[i,])) }

allNames <- grep("allCL", ls(), value=TRUE)

for (cop in allNames) {  
	obj <- get(cop)
	names(obj) <- row.names(obj)
	assign(cop, obj, envir=globalenv())
}

for (name in allNames) { makeState(cptree, get(name),(paste0("anc",name))) }



