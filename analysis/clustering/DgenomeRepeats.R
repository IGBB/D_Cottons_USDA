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


# glm/anova for cluster abundance

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


######## ancestral state reconstruction of "very" significant clusters ######## 

cptree <- read.nexus("concat.raxml.tre")
cptree <- drop.tip(cptree, "'D10_8.o'")

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


               
                             
                             
########### relative aging of transposable elements ###########


