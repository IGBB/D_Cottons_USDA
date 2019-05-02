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
library(dplyr)
library(phytools)
library(matrixStats)

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
# [1] matrixStats_0.52.2 phytools_0.6-20    maps_3.2.0         dplyr_0.7.4        geiger_2.0.6      
# [6] testit_0.7         gridExtra_2.3      reshape2_1.4.2     factoextra_1.0.5   ggrepel_0.7.0     
# [11] scales_0.5.0       ggplot2_2.2.1      geomorph_3.0.5     ape_4.1            rgl_0.98.1       factoextra_1.0.4     ggrepel_0.6.5        plotrix_3.6-4        BiocInstaller_1.22.3 scales_0.4.1         ggplot2_2.2.1       
# 
# loaded via a namespace (and not attached):
#  [1] Rcpp_0.12.10     digest_0.6.12    ggpubr_0.1.2     grid_3.3.3       plyr_1.8.4       gtable_0.2.0     magrittr_1.5     stringi_1.1.5    lazyeval_0.2.0   labeling_0.3     tools_3.3.3     
# [12] stringr_1.2.0    munsell_0.4.3    colorspace_1.3-2 tibble_1.3.0    
# 
##############################################

data <- read.table("combineLibrary.comparativeAnalysis.txt", header = T, sep="\t")
data$size <- as.numeric(rowSums(data[,-1]))
data$percent <- cumsum(data$size)/sum(data$size)
#data <- data[c(300:nrow(data)),]
ggplot(data, aes(x=cluster, y=percent)) + geom_line(size=1) + geom_vline(xintercept=350, color='yellow3', size=1) + scale_x_log10(labels=comma) + scale_y_log10() + geom_vline(xintercept=0, color="grey")+ geom_hline(yintercept=0, color="grey")

Ddata <- read.table("combineLibrary.comparativeAnalysis.txt", header = T, sep="\t")
Ddata <- Ddata[,c(1, 9:(ncol(Ddata)-1))]
Ddata$size <- as.numeric(rowSums(Ddata[,-1]))
Ddata$percent <- cumsum(Ddata$size)/sum(Ddata$size)
#Ddata <- Ddata[c(300:nrow(Ddata)),]
ggplot(Ddata, aes(x=cluster, y=percent)) + geom_line(size=1) + geom_vline(xintercept=350, color='yellow3', size=1) + scale_x_log10(labels=comma) + scale_y_log10() + geom_vline(xintercept=0, color="grey")+ geom_hline(yintercept=0, color="grey")

annots <- read.table("combineLibrary.RMannotations.txt", header = T, sep="\t")
data <- data[c(1:385),]

annot_clust <- cbind(annots, data)
row.names(annot_clust) <- annot_clust$Cluster
annot_clust <- annot_clust[,-c(1,3,39:40)]
names(annot_clust) <- gsub("_{2,}", "", names(annot_clust))
names(annot_clust) <- gsub("_$", "", names(annot_clust))
names(annot_clust) <- gsub("_", ".", names(annot_clust))
names(annot_clust) <- gsub("D3[.]", "D3", names(annot_clust))
names(annot_clust) <- gsub("D2[.]", "D2-", names(annot_clust))


### PCoA of scaled data ###
## with African cottons

Aford_table <- annot_clust[c(1:350),-1]
AfMb_table <- t(Aford_table*0.0085)

AfGStable <- read.table("genome_sizes_Af.txt", header=T, sep="\t")
row.names(AfGStable) <- AfGStable$Species
AfGStable <- AfGStable[order(match(row.names(AfGStable), row.names(AfMb_table))),]
AfGStable$Species <- NULL

Afperc_table <- sweep(AfMb_table, 1, AfGStable$Size, "/")

Afscaledata <- scale(Afperc_table, center = T)
Afd <- dist(Afscaledata, method = "euclidean")

Afcmdfit <- cmdscale(Afd,eig=TRUE, k=2) # k is the number of dim
x <- Afcmdfit$points[,1]
y <- Afcmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(x, y, labels = row.names(AfMb_table), cex=.7)
points(x, y, pch=19) 

Afcmdpoints <- as.data.frame(Afcmdfit$points)
Afcmdpoints$species <- rownames(Afcmdpoints)
Afcmdpoints$species <- gsub("[.].*", "", Afcmdpoints$species)
Afcmdpoints$subgenus <-Afcmdpoints$species
Afcmdpoints$subgenus <- gsub("D4|D7|D9|D11", "Erioxylum", Afcmdpoints$subgenus)
Afcmdpoints$subgenus <- gsub("D2-1|D2-2|D10", "Caducibracteata", Afcmdpoints$subgenus)
Afcmdpoints$subgenus <- gsub("D3D|D3K", "Integrifolia", Afcmdpoints$subgenus)
Afcmdpoints$subgenus <- gsub("D5", "Austroamericana", Afcmdpoints$subgenus)
Afcmdpoints$subgenus <- gsub("D6", "Selera", Afcmdpoints$subgenus)
Afcmdpoints$subgenus <- gsub("D1|D8", "Houzingenia", Afcmdpoints$subgenus)
Afcmdpoints$subgenus <- gsub("A1|A2", "Gossypium", Afcmdpoints$subgenus)
Afcmdpoints$subgenus <- gsub("F1", "Longiloba", Afcmdpoints$subgenus)
Afcmdpoints$species <- rownames(Afcmdpoints)

# should add R hue here to differentiate colors
png("cotton.GS.Af.ordination.png", 5000, 5000, pointsize=12, res=600)
ggplot(Afcmdpoints, aes(x=V1, y=V2, color=subgenus)) + geom_point(size=1) + xlab("PCoA component 1") + ylab("PCoA component 2") + geom_text_repel(aes(label=species))+stat_ellipse()
dev.off()


### New World cottons only

ord_table <- annot_clust[c(1:385),c(9:(ncol(annot_clust)-1))]


# 0.0085 multiplier represents # reads (x) * 95nt/read * 1 kb/1000nt * 1Mb/1000kb * 100% = # reads * 0.0095 = # Mb in entire genome for that cluster 
Mb_table <- t(ord_table*0.0085)

GStable <- read.table("genome_sizes_combine.txt", header=T, sep="\t")
row.names(GStable) <- GStable$Species
GStable <- GStable[order(match(row.names(GStable), row.names(Mb_table))),]
GStable$Species <- NULL

perc_table <- sweep(Mb_table, 1, GStable$Size, "/")

scaledata <- scale(perc_table, center = T)
d <- dist(scaledata, method = "euclidean")

cmdfit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
x <- cmdfit$points[,1]
y <- cmdfit$points[,2]
plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", main="Metric MDS",	type="n")
text(x, y, labels = row.names(Mb_table), cex=.7)
points(x, y, pch=19) 

cmdpoints <- as.data.frame(cmdfit$points)
cmdpoints$species <- rownames(cmdpoints)
cmdpoints$species <- gsub("[.].*", "", cmdpoints$species)
cmdpoints$subgenus <-cmdpoints$species
cmdpoints$subgenus <- gsub("D4|D7|D9|D11", "Erioxylum", cmdpoints$subgenus)
cmdpoints$subgenus <- gsub("D2-1|D2-2|D10", "Caducibracteata", cmdpoints$subgenus)
cmdpoints$subgenus <- gsub("D3D|D3K", "Integrifolia", cmdpoints$subgenus)
cmdpoints$subgenus <- gsub("D5", "Austroamericana", cmdpoints$subgenus)
cmdpoints$subgenus <- gsub("D6", "Selera", cmdpoints$subgenus)
cmdpoints$subgenus <- gsub("D1|D8", "Houzingenia", cmdpoints$subgenus)
cmdpoints$species <- rownames(cmdpoints)

# should add R hue here to differentiate colors
png("cotton.GS.ordination.png", 5000, 5000, pointsize=12, res=600)
ggplot(cmdpoints, aes(x=V1, y=V2, color=subgenus)) + geom_point(size=1) + xlab("PCoA component 1") + ylab("PCoA component 2") + geom_text_repel(aes(label=species))+stat_ellipse()
dev.off()

### PCoA of log-transformed data ###
# evaluates how close the overall repetitive profiles are to one another

ord_table <- ord_table[!(rowSums(ord_table==0)),]
#313 rows left

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
logcmdpoints$species <- rownames(logcmdpoints)
logcmdpoints$species <- gsub("[.].*", "", logcmdpoints $species)
logcmdpoints$subgenus <-logcmdpoints$species
logcmdpoints$subgenus <- gsub("D4|D7|D9|D11", "Erioxylum", cmdpoints$subgenus)
logcmdpoints$subgenus <- gsub("D2-1|D2-2|D10", "Caducibracteata", cmdpoints$subgenus)
logcmdpoints$subgenus <- gsub("D3D|D3K", "Integrifolia", cmdpoints$subgenus)
logcmdpoints$subgenus <- gsub("D5", "Austroamericana", cmdpoints$subgenus)
logcmdpoints$subgenus <- gsub("D6", "Selera", cmdpoints$subgenus)
logcmdpoints$subgenus <- gsub("D1|D8", "Houzingenia", cmdpoints$subgenus)
logcmdpoints$species <- rownames(logcmdpoints)

png("cotton.raw.ordination.log.png", 5000, 5000, pointsize=12, res=600)
ggplot(logcmdpoints, aes(x=V1, y=V2, color=subgenus)) + geom_point(size=2) + xlab("PCoA component 1") + ylab("PCoA component 2") + geom_text_repel(aes(label=species))+stat_ellipse()
dev.off()


### PCA to get the variation ###

cluster.pca <- prcomp(logdata, scale = TRUE)
cluster.eig <- get_eigenvalue(cluster.pca)

ind.coord <- cluster.pca$x

subgen <- as.factor(logcmdpoints$subgenus)

png("cotton_outgroup.PCA.raw.png", 5000, 5000, pointsize=12, res=600)
fviz_pca_ind(cluster.pca, habillage=subgen, pointsize =2, invisible="quali", repel=TRUE, labelsize=3) + theme_minimal() + labs(title = "PCA of log counts") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text(face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5), legend.position="none") +theme_set(theme_grey(base_size=12))
dev.off()

### Procrustea ANOVA
### which species differ from one another in repeats
### https://groups.google.com/forum/#!topic/geomorph-r-package/8_B2thhxU_o

# Adams, D.C., and E. Otarola-Castillo. 2013. geomorph: an R package for the collection and analysis of geometric morphometric shape data. Methods in Ecology and Evolution. 4:393-399.
# Adams, D.C., M. L. Collyer, A. Kaliontzopoulou, and E. Sherratt. 2016 geomorph: Software for geometric morphometric analyses. R package version 3.0.2 http://cran.r-project.org/web/packages/geomorph/index.html.

# procD.lm with data as raw count data
mycountdata <- as.matrix(t(ord_table))
subgen.diff <- advanced.procD.lm(mycountdata~subgen, ~1, groups=~subgen)

#Effect sizes (Z)
#                Austroamericana Caducibracteata   Erioxylum Houzingenia Integrifolia     Selera
#Austroamericana       0.0000000       2.2552632  1.61636436  0.75603405   1.46578041 1.28259878
#Caducibracteata       2.2552632       0.0000000  0.59407635  1.72099640   0.63515443 1.72180713
#Erioxylum             1.6163644       0.5940764  0.00000000  1.29571546  -0.04818612 1.49098921
#Houzingenia           0.7560341       1.7209964  1.29571546  0.00000000   1.23252958 0.01477937
#Integrifolia          1.4657804       0.6351544 -0.04818612  1.23252958   0.00000000 1.83323339
#Selera                1.2825988       1.7218071  1.49098921  0.01477937   1.83323339 0.00000000
#
#P-values
#                Austroamericana Caducibracteata Erioxylum Houzingenia Integrifolia Selera
#Austroamericana           1.000           0.036     0.073       0.206        0.088  0.106
#Caducibracteata           0.036           1.000     0.236       0.067        0.242  0.055
#Erioxylum                 0.073           0.236     1.000       0.104        0.457  0.085
#Houzingenia               0.206           0.067     0.104       1.000        0.118  0.421
#Integrifolia              0.088           0.242     0.457       0.118        1.000  0.046
#Selera                    0.106           0.055     0.085       0.421        0.046  1.000

Afmycountdata <- as.matrix(t(Aford_table))
Afsubgen <- as.factor(Afcmdpoints$subgenus)
Afsubgen.diff <- advanced.procD.lm(Afmycountdata~Afsubgen, ~1, groups=~Afsubgen)

#Effect sizes (Z)
                Austroamericana Caducibracteata  Erioxylum Gossypium Houzingenia Integrifolia
#Austroamericana       0.0000000      -0.8960916 -0.9907966  6.396991  -1.1563209   -1.0412176
#Caducibracteata      -0.8960916       0.0000000 -1.1238477  5.631722  -0.9662491   -1.1049550
#Erioxylum            -0.9907966      -1.1238477  0.0000000  6.228533  -1.1153310   -1.2293992
#Gossypium             6.3969909       5.6317218  6.2285326  0.000000   5.3269414    4.8565614
#Houzingenia          -1.1563209      -0.9662491 -1.1153310  5.326941   0.0000000   -1.0078462
#Integrifolia         -1.0412176      -1.1049550 -1.2293992  4.856561  -1.0078462    0.0000000
#Longiloba             0.7782894       0.6965394  0.7129289  1.156752   0.6288159    0.6002283
#Selera               -1.0797504      -0.9771330 -1.0628491  4.079966  -1.1651657   -0.8649818
#                Longiloba     Selera
#Austroamericana 0.7782894 -1.0797504
#Caducibracteata 0.6965394 -0.9771330
#Erioxylum       0.7129289 -1.0628491
#Gossypium       1.1567521  4.0799659
#Houzingenia     0.6288159 -1.1651657
#Integrifolia    0.6002283 -0.8649818
#Longiloba       0.0000000  0.4800840
#Selera          0.4800840  0.0000000

#P-values
#                Austroamericana Caducibracteata Erioxylum Gossypium Houzingenia Integrifolia
#Austroamericana           1.000           0.831     0.864     0.001       0.951        0.889
#Caducibracteata           0.831           1.000     0.942     0.001       0.832        0.917
#Erioxylum                 0.864           0.942     1.000     0.001       0.929        0.948
#Gossypium                 0.001           0.001     0.001     1.000       0.001        0.001
#Houzingenia               0.951           0.832     0.929     0.001       1.000        0.847
#Integrifolia              0.889           0.917     0.948     0.001       0.847        1.000
#Longiloba                 0.202           0.229     0.203     0.154       0.221        0.254
#Selera                    0.909           0.839     0.889     0.006       0.897        0.747
#                Longiloba Selera
#Austroamericana     0.202  0.909
#Caducibracteata     0.229  0.839
#Erioxylum           0.203  0.889
#Gossypium           0.154  0.006
#Houzingenia         0.221  0.897
#Integrifolia        0.254  0.747
#Longiloba           1.000  0.294
#Selera              0.294  1.000
 

########### characterize composition ###########

# 8.5 multiplier represents # reads (x) * 85nt/read * 1 kb/1000nt * 100% = # reads * 8.5 = # Kb in entire genome for that class 
Kbamount <- data.frame(annot_clust[1], apply(annot_clust[2:ncol(annot_clust)], 2, function (x) x*8.5))
KBsum <- aggregate(. ~Lineage, data=Kbamount, FUN=sum)
names(KBsum) <- gsub("D2\\.", "D2-", names(KBsum))

for (species in unique(gsub("[.].*", "", names(KBsum[-1]))) ) {
    columns <- grep(species, names(KBsum))
    speciesSub <- subset(KBsum, select=as.numeric(columns))
    KBsum[,paste0(species,"mean")] <- rowMeans(speciesSub)
    KBsum[,paste0(species,"min")] <- apply(speciesSub,1,min)
    KBsum[,paste0(species,"max")] <- apply(speciesSub,1,max)
}

KBm <- melt(subset(KBsum, select=grep("mean|Lineage", names(KBsum))))
KBmin <- melt(subset(KBsum, select=grep("min", names(KBsum))))
KBmax <- melt(subset(KBsum, select=grep("max", names(KBsum))))

KBm$min <- KBmin$value
KBm$max <- KBmax$value

limits <- aes(ymax=KBm$max, ymin=KBm$min)

dodge <- position_dodge(width=0.9)

png("Figure_TE.amounts.png", 7500, 5000, pointsize=12, res=600)
ggplot(KBm, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + geom_errorbar(limits, position = dodge) + labs(title = "Aggregate amounts in each species", x="Broad element category", y="Aggregate amount (mean) in kilobases") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text( face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5))+theme_set(theme_grey(base_size=12))
dev.off()

KBsumD <- subset(KBsum, select=grep("D|Lineage", names(KBsum)))
KBmD <- melt(subset(KBsumD, select=grep("mean|Lineage", names(KBsumD))))

KBminD <- melt(subset(KBsumD, select=grep("min", names(KBsumD))))
KBmaxD <- melt(subset(KBsumD, select=grep("max", names(KBsumD))))

KBmD$min <- KBminD$value
KBmD$max <- KBmaxD$value

limits <- aes(ymax=KBmD$max, ymin=KBmD$min)

png("Figure_TE.amounts.Donly.png", 7500, 5000, pointsize=12, res=600)
ggplot(KBmD, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + geom_errorbar(limits, position = dodge) + labs(title = "Aggregate amounts in each species", x="Broad element category", y="Aggregate amount (mean) in kilobases") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text( face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5))+theme_set(theme_grey(base_size=12))
dev.off()





ttable <- 100*(t(Afperc_table))
lin <- Kbamount[,1:2]
lin$A1.155 <- NULL

GSamount <- merge(lin,ttable, by="row.names")
row.names(GSamount) <- GSamount$Row.names
GSamount$Row.names <- NULL

GSsum <- aggregate(. ~Lineage, data=GSamount, FUN=sum)

for (species in unique(gsub("[.].*", "", names(GSsum[-1]))) ) {
    GScolumns <- grep(species, names(GSsum))
    GSspeciesSub <- subset(GSsum, select=as.numeric(GScolumns))
    GSsum[,paste0(species,"mean")] <- rowMeans(GSspeciesSub)
    GSsum[,paste0(species,"min")] <- apply(GSspeciesSub,1,min)
    GSsum[,paste0(species,"max")] <- apply(GSspeciesSub,1,max)
}

GSm <- melt(subset(GSsum, select=grep("mean|Lineage", names(GSsum))))
GSmin <- melt(subset(GSsum, select=grep("min", names(GSsum))))
GSmax <- melt(subset(GSsum, select=grep("max", names(GSsum))))

GSm$min <- GSmin$value
GSm$max <- GSmax$value

limits <- aes(ymax=GSm$max, ymin=GSm$min)


png("Figure_TE.amounts.GS.png", 7500, 5000, pointsize=12, res=600)

ggplot(GSm, aes(x=Lineage, y=value, fill = variable)) + geom_bar(stat = "identity",position = dodge) + geom_errorbar(limits, position = dodge) + labs(title = "Percent Genome Size", x="Broad element category", y="Percent Genome Size") + theme(axis.text = element_text(size = rel(1.5)), plot.margin=margin(2,2,2,2,"cm"), plot.title=element_text( face="bold", hjust=0.5), axis.title.x = element_text(face="bold", hjust=0.5), axis.title.y = element_text(face="bold", vjust=0.5))+theme_set(theme_grey(base_size=12))

dev.off()


########################TE Amounts ########################



TEamounts <- as.data.frame(cbind(TEtotal=c(rowSums(Mb_table)), GS=c(as.matrix(GStable[,1])), TEperc=c((rowSums(Mb_table)/as.matrix(GStable[,1])))))

gClust <- rownames(annot_clust)[annot_clust$Lineage=="LTR/Gypsy"]
gTable <- Mb_table[,gClust]
gTEamounts <- as.data.frame(cbind(gTable=c(rowSums(gTable)), GS=c(as.matrix(GStable[,1])), TEperc=c((rowSums(gTable)/as.matrix(GStable[,1])))))



# TEamounts[which.min(TEamounts$TEperc),]
#      TEtotal  GS  TEperc
# D5.53  311.508 880 0.3539864

# TEamounts[which.max(TEamounts$TEperc),]
#       TEtotal  GS    TEperc
# D5.6 509.49 880 0.5789659

# TEamounts[which.max(TEamounts$TEperc[!row.names(TEamounts)=="D5.6"]),]
#       TEtotal  GS    TEperc
# D6.5 351.798 841 0.4183092


# gTEamounts[which.min(gTEamounts$TEperc),]
#       gTEtotal  GS    TEperc
# D5.53  225.097 880 0.255792

# gTEamounts[which.max(gTEamounts$TEperc),]
#      gTEtotal  GS    TEperc
# D5.6  363.5705 880 0.4131483

# gTEamounts[which.max(gTEamounts$TEperc[!row.names(gTEamounts)=="D5.6"]),]
#      gTEtotal  GS    TEperc
# D6.5  267.0955 841 0.3175927

TEav <- TEamounts
TEav$species <- gsub("\\..*", "", row.names(TEav))

av <- aggregate(. ~species, data=TEav, FUN=mean)
av[which.min(av$TEperc),]
# species  TEtotal  GS    TEperc
#    D2-2  358.989  910 0.3944934

av[which.max(av$TEperc),]
# species  TEtotal  GS    TEperc
#    D2-1 401.1575 856 0.4686419


gTEav <- gTEamounts
gTEav$species <- gsub("\\..*", "", row.names(gTEav))

gav <- aggregate(. ~species, data=gTEav, FUN=mean)
gav[which.min(gav$TEperc),]
# species  TEtotal  GS    TEperc
#     D10  D10 266.05 910 0.2923626

gav[which.max(gav$TEperc),]
# species  TEtotal  GS    TEperc
#      D6 288.7833 841 0.3433808


rowMeans(Mbsum[7,c(2:8,36)])/1000
# 861.0524  ## Mb of gypsy in 385 clusters

rowMeans(Mbsum[7,c(9:35)])/1000
# 276.8592  ## Mb of gypsy in 385 clusters


rowMeans(Mbsum[4,c(2:8,36)])/1000
# 1.558688  ## Mb of MuDR in 385 clusters

rowMeans(Mbsum[4,c(9:35)])/1000
# 4.390407  ## Mb of MuDr in 385 clusters


rowMeans(Mbsum[6,c(2:8)])/1000
# 41.344   ## Mb of copia in 385 clusters

rowMeans(Mbsum[6,c(9:35)])/1000
# 37.4192  ## Mb of copia in 385 clusters

(Mbsum[6,36])/1000
# 39.3975  ## Mb of copia in 385 clusters


##########################################################

getGLM <- function (df, factors)
{
    data.in <- t(df)

    P.scores <- rep(NA, ncol(data.in))
    for (i in 1:ncol(data.in)){
        P.scores[i] <- anova(lm(data.in[,i]~factors))$`Pr(>F)`[1]  
    }

    df$p.value <- P.scores

    df$p.valueBH <- p.adjust(df$p.value, method="BH")

    # get rid of NA
    df$p.valueBH[is.na(df$p.valueBH)] <- 99
    return(df)
}


# glm/anova for cluster abundance

GLMclust <- subset(annot_clust, select=grep("A|D", names(annot_clust)))
subsections <- as.factor(gsub("[1234567890].*", "", names(GLMclust)))

GLMD <- subset(annot_clust, select=grep("D", names(annot_clust)))
subD <- as.factor(gsub("D3.*", "In", names(GLMD)))
subD <- as.factor(gsub("D5.*", "Au",subD))
subD <- as.factor(gsub("D6.*", "Se", subD))
subD <- as.factor(gsub("D11.*", "Er", subD))
subD <- as.factor(gsub("D[479].*", "Er", subD))
subD <- as.factor(gsub("D[18].*", "Ho", subD))
subD <- as.factor(gsub("D[[21][-0].*", "Ca", subD))

GLMDn6 <- GLMD[,-c(18)]
subDn6 <- subD[-18]

GLMclust <- getGLM(GLMclust,subsections)
GLMD <- getGLM(GLMD,subD)
GLMDn6 <- getGLM(GLMDn6,subDn6)

### how many clusters differ between A and D
length(row.names(GLMclust[(GLMclust$p.valueBH<0.001),]))
# 234 out of 385 clusters are significantly different among these

### how many clusters differ among Ds
length(row.names(GLMD[(GLMD$p.valueBH<0.05),]))
# 14 out of 385 clusters are significantly different between these

### how many clusters differ among Ds
length(row.names(GLMDn6[(GLMDn6$p.valueBH<0.05),]))
# 19 out of 385 clusters are significantly different between these


# glm/anova for D6 cluster abundance

GLMDnoD6 <- subset(annot_clust, select=grep("D", names(annot_clust)))
subD <- as.factor(gsub("D3.*", "In", names(GLMD)))
subD <- as.factor(gsub("D5.*", "Au",subD))
subD <- as.factor(gsub("D6.*", "Se", subD))
subD <- as.factor(gsub("D11.*", "Er", subD))
subD <- as.factor(gsub("D[479].*", "Er", subD))
subD <- as.factor(gsub("D[18].*", "Ho", subD))
subD <- as.factor(gsub("D[[21][-0].*", "Ca", subD))

GLMDnoD6 <- GLMDnoD6[,-c(18)]
subDnoD6 <- gsub("Ho|Au|Ca|Er|In", "other", subDn6)
GLMDnoD6 <- getGLM(GLMDnoD6,subDnoD6)

length(row.names(GLMDnoD6[(GLMDnoD6$p.valueBH<0.001),]))
# 0 out of 385 clusters are significantly different among these




GLMAD6 <- GLMclust %>% select(A1.155, A1.97, A2.1011, A2.34, A2.4, A2.44, A2.JCVI, D6.5, D6.7)
subAD6 <- subDnoD6[c(1:7, 20:21)]
GLMAD6 <- getGLM(GLMAD6,subAD6)

length(row.names(GLMAD6[(GLMAD6$p.valueBH<0.001),]))
# 78 out of 385 clusters are significantly different among these


GLMAD8 <- GLMclust %>% select(A1.155, A1.97, A2.1011, A2.34, A2.4, A2.44, A2.JCVI, D8.8, D8.9)
subAD8 <- subDnoD6[c(1:7, 20:21)]
GLMAD8 <- getGLM(GLMAD8,subAD8)

length(row.names(GLMAD8[(GLMAD8$p.valueBH<0.001),]))
# 79 out of 385 clusters are significantly different among these


GLMAD1 <- GLMclust %>% select(A1.155, A1.97, A2.1011, A2.34, A2.4, A2.44, A2.JCVI, D1.2, D1.35)
subAD1 <- subDnoD6[c(1:7, 20:21)]
GLMAD1 <- getGLM(GLMAD1,subAD1)

length(row.names(GLMAD1[(GLMAD1$p.valueBH<0.001),]))
# 68 out of 385 clusters are significantly different among these



### let's grab the clusters that differentiate the Ds

clustersD <- row.names(GLMD[(GLMD$p.valueBH<0.05),])

sigD <- annot_clust[clustersD,]

table(sigD$Lineage)
#        *           DNA DNA/MULE-MuDR           LTR     LTR/Copia     LTR/Gypsy          rRNA   Unspecified 
#        1             0             0             0             0            13             0             0 

table(annot_clust$Lineage)
#  *           DNA       DNA/hAT DNA/MULE-MuDR           LTR     LTR/Copia     LTR/Gypsy          rRNA   Unspecified 
# 29             1             1             3            40            43           262             1             5 



nrow(annot_clust[annot_clust$Lineage=="LTR/Gypsy",])/nrow(annot_clust)
# [1] 0.6805195
nrow(annot_clust[annot_clust$Lineage=="LTR/Copia",])/nrow(annot_clust)
# [1] 0.1116883
nrow(annot_clust[annot_clust$Lineage=="LTR",])/nrow(annot_clust)
# [1] 0.1038961

readsD <- GLMD[(GLMD$p.valueBH<0.05),]



######## ancestral state reconstruction of 3 (each up and down) randomly sampled "very" significant clusters ######## 

ancClusters <- row.names(readsD)
# [1] "CL0030" "CL0053" "CL0074" "CL0086" "CL0142" "CL0154" "CL0167" "CL0193" "CL0201" "CL0301"
# [11] "CL0339" "CL0341" "CL0344"  "CL0357"

cptree <- read.nexus("concat.raxml.tre")
cptree <- drop.tip(cptree, "'D10_8.o'")
cptree <- drop.tip(cptree, "D4_12C")

plot.phylo(cptree, cex =1, label.offset=0.0005, align.tip.label=TRUE)
nodelabels(cptree$node, adj=c(1.1,-0.2), frame="none", cex=1)

ancTable <- Kbamount[ancClusters,-c(1:8)]
ancTable <- ancTable[, -which(names(ancTable) %in% c("D4.12"))]

### match names####
names(ancTable) <- gsub("\\.", "_", names(ancTable))
names(ancTable)[names(ancTable) == 'D5ref'] <- 'D5_ref'
names(ancTable)[names(ancTable) == 'D3D_2'] <- 'D3D_27'
names(ancTable)[names(ancTable) == 'D3K_5'] <- 'D3K_57'
names(ancTable)[names(ancTable) == 'D4_12'] <- 'D4_12C'
### end match names ###



for (i in c(1:length(ancClusters))) { assign(paste0("anc",ancClusters[i]), t(subset(ancTable, row.names(ancTable) %in% ancClusters[i]))) }


varNames <- grep("anc", ls(), value=TRUE)
rmAn <- c("ancClusters", "ancTable")
varNames <- varNames[!varNames %in% rmAn]

for (anc in varNames) {  
	obj <- get(anc)
	names(obj) <- row.names(obj)
	assign(anc, obj, envir=globalenv())
}


name.check(cptree,ancCL0030) # check one table to make sure the names match
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

for (name in varNames) { checkModel(cptree, get(name),outdf=paste0("mdf",name)) }

mdfNames <- grep("mdfancCL", ls(), value=TRUE)
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

# if lambda factor is so close to 1 then ~ BM anyway

# example rescales 
# delta.anc360 <- fitContinuous(cptree, anc360 , model="delta")
# delta.tree360 <- rescale(cptree, "delta", delta.anc360$opt$delta)
# kappa.anc96 <- fitContinuous(cptree, anc96 , model="kappa")
# kappa.tree96 <- rescale(cptree, "kappa", kappa.anc96$opt$kappa)
# lambda.anc20 <- fitContinuous(cptree, anc20, model="lambda")
# lambda.tree20 <- rescale(cptree, "lambda", lambda.anc20$opt$lambda)

bestModel
#        cluster bestModel         modelAIC            BMAIC     diff          exp
#1  mdfancCL0030        DE 313.074065163735 333.609235198952 20.53517 2.878429e+04
#14 mdfancCL0357        DE 366.202739841951   368.6111710846  2.408431 3.334143e+00
#2  mdfancCL0053        KA 442.451746384012 465.495956417128 23.04421 1.009222e+05
#4  mdfancCL0086        KA 437.966272213003 480.297724909138 42.33145 1.556532e+09
#5  mdfancCL0142        KA 421.127178953575 432.387809694809 11.26063 2.787500e+02
#8  mdfancCL0193        KA 469.086576007268 488.399835414019 19.31326 1.562503e+04
#10 mdfancCL0301        KA  409.65753048415 490.897944607258 81.24041 4.376549e+17
#11 mdfancCL0339        KA 368.174380772029 424.220877737801 56.04650 1.480274e+12
#13 mdfancCL0344        KA 394.680865331074 458.983240721868 64.30238 9.185090e+13
#3  mdfancCL0074        LA 443.533242914413  494.19920247773 50.66596 1.004554e+11
#6  mdfancCL0154        LA 402.002759639846 429.946864955597 27.94411 1.169460e+06
#7  mdfancCL0167        LA 392.127015985954 438.159865518199 46.03285 9.906181e+09
#9  mdfancCL0201        LA 405.177761711182 433.137498607341 27.95974 1.178636e+06
#12 mdfancCL0341        LA 312.459543043772 378.510010677098 66.05047 2.201288e+14


delta.ancCL0030 <- fitContinuous(cptree, ancCL0030 , model="delta")
delta.treeCL0030 <- rescale(cptree, "delta", delta.ancCL0030$opt$delta)

delta.ancCL0357 <- fitContinuous(cptree, ancCL0357 , model="delta")
delta.treeCL0357 <- rescale(cptree, "delta", delta.ancCL0357$opt$delta)

for (i in c("0053", "0086", "0142", "0193", "0301", "0339", "0344")){
    x <- get(paste0("ancCL",i))
    kanc <- fitContinuous(cptree, x , model="kappa")
    ktree <- rescale(cptree, "kappa", kanc$opt$kappa)
    assign(paste0("kappa.ancCL",i), kanc, envir=globalenv())
    assign(paste0("kappa.treeCL",i), ktree, envir=globalenv())
}

for (i in c("0074", "0154", "0167", "0201", "0341")){
    x <- get(paste0("ancCL",i))
    lanc <- fitContinuous(cptree, x , model="lambda")
    ltree <- rescale(cptree, "lambda", lanc$opt$lambda)
    assign(paste0("lambda.ancCL",i), lanc, envir=globalenv())
    assign(paste0("lambda.treeCL",i), ltree, envir=globalenv())
}


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
    GSfit <- fastAnc(tree, df, vars=TRUE, CI=TRUE)
    gfit <- round(GSfit$ace)
    assign(paste0("G",name), GSgradient, envir=globalenv())
    assign(paste0("fit",name),gfit,envir=globalenv())
}


makeState(delta.treeCL0030, ancCL0030, name="ancCL0030")
makeState(delta.treeCL0357, ancCL0357, name="ancCL0357")

for (i in c("0053", "0086", "0142", "0193", "0301", "0339", "0344")){
    makeState(get(paste0("kappa.treeCL",i)), get(paste0("ancCL",i)),(paste0("ancCL",i))) }

for (i in c("0074", "0154", "0167", "0201", "0341")){
    makeState(get(paste0("lambda.treeCL",i)), get(paste0("ancCL",i)),(paste0("ancCL",i))) }

nodeSize <- ls(pattern="fitanc")

vals <- vector(mode="numeric", length=0) # compare tips to reconstructed basal D node
Fvals <- vector(mode="numeric", length=0)  # compare tips to reconstructed basal node
Fchange <- vector(mode="numeric", length=0)

for (name in row.names(ancCL0053)) { 
    newvec <- vector(mode="numeric", length=0)
    assign(paste0("vec",name),newvec,envir=globalenv()) }


for (i in varNames) { 
    
    temp <- get(i)
    sub <- get(paste0("fit",i))

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

boxplot.matrix <- cbind(vecD1_2, vecD1_35, vecD2_1_6, vecD2_2, vecD3D_27, vecD3K_56, vecD3K_57, vecD4_185, vecD5_2, vecD5_31, vecD5_4, vecD5_53, vecD5_6, vecD5_8, vecD5_ref, vecD6_5, vecD6_7, vecD7_157, vecD7_4, vecD8_8, vecD8_9, vecD9_4, vecD10_3, vecD10_7, vecD10_8, vecD11_1)
bpm <- melt(boxplot.matrix)
bpm$sub <- bpm$Var2
bpm$sub <- gsub("vecD6_.*","Selera",bpm$sub)
bpm$sub <- gsub("vecD5_.*", "Austroamericana",bpm$sub)
bpm$sub <- gsub("vecD4_.*|vecD7_.*|vecD9_.*|vecD11_.*","Erioxylum",bpm$sub)
bpm$sub <- gsub("vecD3.*","Integrifolia",bpm$sub)
bpm$sub <- gsub("vecD2.*|vecD10_.*", "Caducibracteata",bpm$sub)
bpm$sub <- gsub("vecD1_.*|vecD8_.*", "Houzingenia",bpm$sub)

png("Figure.magnitude.gainloss.diffCount.png", 15000, 7320, pointsize=12, res=600)
ggplot(bpm, aes(x=Var2, y=value, fill=sub)) + geom_violin()+scale_y_log10()+geom_boxplot(width=0.1, fill="black")+coord_flip()
dev.off()

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

write.table(plusminus, file="gainloss.tbl", quote=F, sep="\t", row.names=TRUE)

length(which(Fvals<0))
#[1] 111

length(which(Fvals>0))
#[1] 240

median(Fvals)
#[1] 247.5


gainloss <- c(sum(vecD1_2), sum(vecD1_35), sum(vecD2_1_6), sum(vecD2_2), sum(vecD3D_27), sum(vecD3K_56), sum(vecD3K_57), sum(vecD4_185), sum(vecD5_2), sum(vecD5_31), sum(vecD5_4), sum(vecD5_53), sum(vecD5_6), sum(vecD5_8), sum(vecD5_ref), sum(vecD6_5), sum(vecD6_7), sum(vecD7_157), sum(vecD7_4), sum(vecD8_8), sum(vecD8_9), sum(vecD9_4), sum(vecD10_3), sum(vecD10_7), sum(vecD10_8), sum(vecD11_1))
names(gainloss)<- c("D1_2", "D1_35", "D2_1_6", "D2_2", "D3D_27", "D3K_56", "D3K_57", "D4_185", "D5_2", "D5_31", "D5_4", "D5_53", "D5_6", "D5_8", "D5_ref", "D6_5", "D6_7", "D7_157", "D7_4", "D8_8", "D8_9", "D9_4", "D10_3", "D10_7", "D10_8", "D11_1")

gainloss 
#   D1_2   D1_35  D2_1_6    D2_2  D3D_27  D3K_56  D3K_57  D4_185    D5_2   D5_31    D5_4   D5_53    D5_6    D5_8  D5_ref 
# -369.0   549.0  4739.5  4297.5  7425.5 12075.0 12262.0   923.0 -2902.0 -2681.0 -3208.0 -2995.5  6524.5 -2502.5   285.5 
#   D6_5    D6_7  D7_157    D7_4    D8_8    D8_9    D9_4   D10_3   D10_7   D10_8   D11_1 
# 3838.5  4595.0  4272.0  5530.0  3039.5  1305.5  3039.5  4374.0  4841.5  4204.0   880.5 

gainloss <- melt(gainloss)

png("Figure.magnitude.gainloss.diffCount.png", 15000, 7320, pointsize=12, res=600)
ggplot(gainloss, aes(x=row.names(gainloss),y=value)) + geom_bar(stat="identity")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

png("Figure.number.gainloss.png", 15000, 7320, pointsize=12, res=600)
gldf <- as.data.frame(cbind(vals, Fvals, Fchange))
ggplot(gldf, aes(x=vals)) + geom_density(fill="blue",color="blue")+geom_vline(xintercept=0, size=1, linetype="dashed")
dev.off()                             
         







ls(pattern="Ganc")
#  [1] "GancCL0030" "GancCL0053" "GancCL0074" "GancCL0086" "GancCL0142" "GancCL0154" "GancCL0167" "GancCL0193" "GancCL0201"
# [10] "GancCL0301" "GancCL0339" "GancCL0341" "GancCL0344" "GancCL0357" 




png("Figure_grid.anc.14.png", 10000, 7500, pointsize=12, res=600)
split.screen(figs=c(4,4))

screen(1)

par(fg="transparent", bg="white")
plot(GancCL0030, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0030), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="white")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0030$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(2)
par(fg="transparent")
plot(GancCL0053, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0053), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0053$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")

screen(3)
par(fg="transparent")
plot(GancCL0074, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0074), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0074$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(4)
par(fg="transparent")
plot(GancCL0086, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0086), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0086$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(5)
par(fg="transparent")
plot(GancCL0142, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0142), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0142$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(6)
par(fg="transparent")
plot(GancCL0154, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0154), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0154$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(7)
par(fg="transparent")
plot(GancCL0167, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0167), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0167$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(8)
par(fg="transparent")
plot(GancCL0193, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0193), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0193$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(9)
par(fg="transparent")
plot(GancCL0201, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0201), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0201$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(10)
par(fg="transparent")
plot(GancCL0301, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0301), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0301$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(11)
par(fg="transparent")
plot(GancCL0339, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0339), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0339$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(12)
par(fg="transparent")
plot(GancCL0341, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0341), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0341$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")


screen(13)
par(fg="transparent")
plot(GancCL0344, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0344), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0344$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")

screen(14)
par(fg="transparent")
plot(GancCL0357, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(fitancCL0357), adj=c(-0.3,0.2), frame="none", cex=0.6)
par(fg="black", bg="transparent")
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
xlim<-obj$x.lim
ylim<-obj$y.lim
offset<-0.1*diff(xlim)
xlim[2]<-xlim[2]+offset
cw<-reorder(GancCL0357$tree)
text(rep(max(obj$xx[1:Ntip(cw)])+offset,Ntip(cw)),obj$yy[1:Ntip(cw)], labels=gsub("_"," ",cw$tip.label),font=3,pos=4,cex=0.6,offset=0.02)
for(i in 1:Ntip(cw)) lines(c(obj$xx[i],max(obj$xx[1:Ntip(cw)])+offset), rep(obj$yy[i],2),lty="dotted")

close.screen(all.screens=TRUE)


dev.off()


########### gypsys in D5 ########### 

D5gyp <- annot_clust[annot_clust$Lineage=="LTR/Gypsy", c(22:28)]
nrow(D5gyp)
# 262


D5gypDiff <- data.frame(accession=character(), Max=numeric(), threeSD=numeric(), twoSD=numeric(), oneSD=numeric(), stringsAsFactors=FALSE)

for (i in c(1:7)) {
    x <- D5gyp
    x$otherMean <- rowMeans(x[,-i])
    x$oneSD <- rowSds(as.matrix(x[,-c(i,8)]), na.rm=TRUE)
    x$twoSD <- 2*rowSds(as.matrix(x[,-c(i,8:9)]), na.rm=TRUE)
    x$threeSD <- 3*rowSds(as.matrix(x[,-c(i,8:10)]), na.rm=TRUE)
    D5gypDiff[i,1] <- names(x[i])
    D5gypDiff[i,2] <- round((100*nrow(x[which(x[,i] > rowMaxs(as.matrix(x[,-i]))), ])/nrow(D5gyp)),2)
    D5gypDiff[i,3] <- round((100*nrow(x[which(x[,i] > x$threeSD), ])/nrow(D5gyp)),2)
    D5gypDiff[i,4] <- round((100*nrow(x[which(x[,i] > x$twoSD), ])/nrow(D5gyp)),2)
    D5gypDiff[i,5] <- round((100*nrow(x[which(x[,i] > x$oneSD), ])/nrow(D5gyp)),2)
}

D5gyp$Mean <- rowMeans(D5gyp)
D5gyp$oneSD2 <- rowSds(as.matrix(D5gyp[,-c(5,8)]), na.rm=TRUE)
D5gyp$twoSD <- 2*rowSds(as.matrix(D5gyp[,-c(5,8:9)]), na.rm=TRUE)
D5gyp$threeSD <- 3*rowSds(as.matrix(D5gyp[,-c(5,8:10)]), na.rm=TRUE)

D5gypDiff["inc",1] <- "D5.6.inclusive"
D5gypDiff["inc",2] <- round((100*nrow(D5gyp[which(D5gyp$D5.6 > rowMaxs(as.matrix(D5gyp[,-c(5,8:10)]))), ])/nrow(D5gyp)),2)
D5gypDiff["inc",3] <- round((100*nrow(D5gyp[which(D5gyp$D5.6 > D5gyp$threeSD), ])/nrow(D5gyp)),2)
D5gypDiff["inc",4] <- round((100*nrow(D5gyp[which(D5gyp$D5.6 > D5gyp$twoSD), ])/nrow(D5gyp)),2)
D5gypDiff["inc",5] <- round((100*nrow(D5gyp[which(D5gyp$D5.6 > D5gyp$oneSD), ])/nrow(D5gyp)),2)

write.table(D5gypDiff, file="D5.gypsy.clusters.tbl", sep="\t", quote=F, row.names=F)

excess <- D5gyp[which(D5gyp$D5.6 > rowMaxs(as.matrix(D5gyp[,-c(5,8:10)]))), c(1:7)]
excess$Max <- round(rowMaxs(as.matrix(excess[,-c(5,8:10)])),0)
excess$diff <- excess$D5.6 - excess$Max
excess$diffMb <- excess$diff*0.0085

nrow(excess)
# [1] 73 
nrow(excess[which(excess$diffMb >=1),])
# [1] 18

rowMeans(excess["CL0078",c(1:7)])*0.0085
#   CL0078 
# 5.226286 



########### copia in D ########### 

copia <- annot_clust[annot_clust$Lineage=="LTR/Copia", ]
nrow(copia)
# 43

copSum <- colSums(copia[,-1])*0.0085  # Make it into Mb
# A1.155   A1.97 A2.1011   A2.34    A2.4   A2.44 A2.JCVI   D10.3   D10.7   D10.8   D11.1    D1.2   D1.35  D2-1.6    D2-2  D3D.27  D3K.56  D3K.57   D4.12  D4.185    D5.2   D5.31 
#41.5395 45.3475 35.2325 39.5505 41.7775 44.6080 41.3525 38.1735 33.8895 35.1390 36.1590 40.4600 34.5015 40.2475 33.3455 36.7115 33.6600 40.1710 39.9925 35.9635 39.4400 34.0425 
#   D5.4   D5.53    D5.6    D5.8  D5.ref    D6.5    D6.7  D7.157    D7.4    D8.8    D8.9    D9.4    F1.1 
#33.6515 29.8180 57.1880 36.3120 33.6600 34.9690 36.2865 34.4335 40.6810 42.6530 38.5985 40.1710 39.3975

mean(copSum[(1:7)])
# [1] 41.344  A average

mean(copSum[(8:(length(copSum)-1))])
# 37.4192 D aberage

# F1 = 39.3975 Mb

sort(copSum)
#  D5.53    D2-2    D5.4  D3K.56  D5.ref   D10.7   D5.31  D7.157   D1.35    D6.5   D10.8 A2.1011  D4.185   D11.1    D6.7    D5.8  D3D.27   D10.3    D8.9    F1.1    D5.2   A2.34 
#29.8180 33.3455 33.6515 33.6600 33.6600 33.8895 34.0425 34.4335 34.5015 34.9690 35.1390 35.2325 35.9635 36.1590 36.2865 36.3120 36.7115 38.1735 38.5985 39.3975 39.4400 39.5505 
#  D4.12  D3K.57    D9.4  D2-1.6    D1.2    D7.4 A2.JCVI  A1.155    A2.4    D8.8   A2.44   A1.97    D5.6 
#39.9925 40.1710 40.1710 40.2475 40.4600 40.6810 41.3525 41.5395 41.7775 42.6530 44.6080 45.3475 57.1880 


######## ancestral state reconstruction of copia clusters ######## 


copiaclusters <- row.names(copia)
copiaclusters 
#  [1] "CL0039" "CL0082" "CL0085" "CL0089" "CL0093" "CL0098" "CL0110" "CL0118" "CL0122"
# [10] "CL0133" "CL0139" "CL0143" "CL0146" "CL0147" "CL0155" "CL0166" "CL0170" "CL0173"
# [19] "CL0186" "CL0194" "CL0197" "CL0214" "CL0215" "CL0221" "CL0224" "CL0244" "CL0254"
# [28] "CL0256" "CL0260" "CL0274" "CL0281" "CL0289" "CL0292" "CL0299" "CL0302" "CL0318"
# [37] "CL0321" "CL0325" "CL0326" "CL0334" "CL0340" "CL0342"


copTable <- Kbamount[copiaclusters ,-c(1:8)]
names(copTable) <- gsub("\\.", "_", names(copTable))
copTable$D4_12 <- NULL



for (i in c(1:length(copiaclusters))) { assign(paste0("cop",copiaclusters[i]), t(subset(copTable, row.names(copTable) %in% copiaclusters[i]))) }


varNames <- grep("copCL", ls(), value=TRUE)
#varNames <- varNames[c(2:16)]

for (cop in varNames) {  
	obj <- get(cop)
	names(obj) <- row.names(obj)
	assign(cop, obj, envir=globalenv())
}


name.check(cptree,copCL0039) # check one table to make sure the names match
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

for (name in row.names(copCL0039)) { 
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

boxplot.matrix <- cbind(vals, Fvals,vecD1_2, vecD1_35, vecD2_1_6, vecD2_2, vecD3D_27, vecD3K_56, vecD3K_57, vecD4_185, vecD5_2, vecD5_31, vecD5_4, vecD5_53, vecD5_6, vecD5_8, vecD5_ref, vecD6_5, vecD6_7, vecD7_157, vecD7_4, vecD8_8, vecD8_9, vecD9_4, vecD10_3, vecD10_7, vecD10_8, vecD11_1, Fchange)


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
#[1] 560

length(which(Fvals>0))
#[1] 571

median(Fvals)
#[1] 1.75


gainloss <- c(sum(vecD1_2), sum(vecD1_35), sum(vecD2_1_6), sum(vecD2_2), sum(vecD3D_27), sum(vecD3K_56), sum(vecD3K_57), sum(vecD4_185), sum(vecD5_2), sum(vecD5_31), sum(vecD5_4), sum(vecD5_53), sum(vecD5_6), sum(vecD5_8), sum(vecD5_ref), sum(vecD6_5), sum(vecD6_7), sum(vecD7_157), sum(vecD7_4), sum(vecD8_8), sum(vecD8_9), sum(vecD9_4), sum(vecD10_3), sum(vecD10_7), sum(vecD10_8), sum(vecD11_1))
names(gainloss)<- c("D1_2", "D1_35", "D2_1_6", "D2_2", "D3D_27", "D3K_56", "D3K_57", "D4_185", "D5_2", "D5_31", "D5_4", "D5_53", "D5_6", "D5_8", "D5_ref", "D6_5", "D6_7", "D7_157", "D7_4", "D8_8", "D8_9", "D9_4", "D10_3", "D10_7", "D10_8", "D11_1")

#   D1_2   D1_35  D2_1_6    D2_2  D3D_27  D3K_56  D3K_57  D4_185    D5_2   D5_31    D5_4 
# 1936.0 -3869.5  1910.5 -4991.5 -1625.5 -4694.0  1834.0 -2373.5  1018.0 -4320.0 -4762.0 
#  D5_53    D5_6    D5_8  D5_ref    D6_5    D6_7  D7_157    D7_4    D8_8    D8_9    D9_4 
#-8612.5 18528.0 -2229.0 -4779.0 -3351.0 -2127.0 -3929.0  2293.0  4256.5   108.5  1757.5 
#  D10_3   D10_7   D10_8   D11_1 
# -180.5 -4490.0 -3249.0 -2246.0 

 
gainloss <- melt(gainloss)

png("Figure.magnitude.gainloss.copia.png", 15000, 7320, pointsize=12, res=600)
ggplot(gainloss, aes(x=row.names(gainloss),y=value)) + geom_bar(stat="identity")+ theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
dev.off()

png("Figure.number.gainloss.copia.png", 15000, 7320, pointsize=12, res=600)
gldf <- as.data.frame(cbind(vals, Fvals, Fchange))
ggplot(gldf, aes(x=vals)) + geom_density(fill="blue",color="blue")+geom_vline(xintercept=0, size=1, linetype="dashed")
dev.off()                             
                             


########### relative aging of transposable elements ###########

# see TE_dating_histogram.pl for table generation

ages <- read.table("TE_dating.txt", header=T, row.names=1, sep="\t")
ages$age <- ages$Category
ages <- ages[c(1:385),]

allages <- read.table("TE_dating.table.txt", header =T, row.names=1, sep="\t")
allages <- allages[c(1:385),]

#Category	age
#1	young
#3	old
#4	young
#5	old
#4*	old

ages$age <- sub("3|5|4\\*","old",ages$age)
ages$age <- sub("1|4","young",ages$age)

Dallages <- as.data.frame(sapply(allages, function(x) sub("3|5|4\\*","old",x)))
Dallages <- as.data.frame(sapply(Dallages, function(x) sub("1|4","young",x)))
names(Dallages) <- gsub("D5_", "D5\\.", names(Dallages))
names(Dallages) <- gsub("_", "", names(Dallages))
names(Dallages) <- gsub("D3\\.", "D3", names(Dallages))
Dallages$D10.8o <- NULL
row.names(Dallages) <- row.names(allages)


table(ages$age)
#  old young 
#  301    84 


ageTable <- data.frame(old=numeric(), young=numeric(), percent=numeric(), stringsAsFactors=FALSE)

for (name in names(Dallages)) {
    x <- as.vector(Dallages[,name])
    x <- x[!is.na(x)]
    x <- x[!x %in% c(2,6)]
    ageTable[name,1] <- table(x)[[1]]
    ageTable[name,2] <- table(x)[[2]]
    ageTable[name,3] <- round(100*ageTable[name,1]/(ageTable[name,1]+ageTable[name,2]), 1)
}

min(ageTable$percent)
#[1] 68.6
max(ageTable$percent)
#[1] 78.6

write.table(ageTable, file="D.ages.tbl", quote=F, sep="\t", row.names=TRUE)

DsmYoung <- Dallages[rowSums(Dallages == "young", na.rm = TRUE)>0,]
DsmYoung$number <- rowSums(DsmYoung == "young", na.rm = TRUE)
table(DsmYoung$number)

# 1  2  3  4  5  6  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 
#45 25 14 14 13 12  8  8  6  1  9  2  1  2  3  4  4  5  9  9  5  8  6  5  5  2 

nrow(DsmYoung)
# 225


manyYoung <- as.data.frame(table(DsmYoung$number))

png("Supp3.B.png", 5000, 3750, pointsize=12, res=600)
ggplot(data=manyYoung, aes(x=Var1, y=Freq, group=1)) + geom_area()
dev.off()



youngAndOld <- c(ageTable$old, ageTable$young)
youngAndOldnames <- c(rep("old", times=27), rep("young", times=27))
spec <- c(row.names(ageTable), row.names(ageTable))
vio <- as.data.frame(cbind(youngAndOldnames, youngAndOld, spec))
names(vio) <- c("name", "num", "species")
vio$num <- as.numeric(as.character(vio$num))
vio$species <- gsub("D2\\.", "D2-", vio$species)
vio$species <- gsub("\\..*", "", vio$species)

png("Supp3.A.png", 5000, 5000, pointsize=12, res=600)
ggplot(vio, aes(factor(name), num, color=species)) + geom_jitter(height=0, width=0.05) + xlab("cluster age") + ylab("number of clusters")
dev.off()










png("FigureS3.hist.number.young.png", 15000, 7320, pointsize=12, res=600)
vioplot(DsmYoung$number, col="light blue")
dev.off()



youngClusters <- as.character(row.names(DsmYoung))

youngTable <- Kbamount[youngClusters ,-c(1:8)]
names(youngTable)<- gsub("\\.", "_", names(youngTable))
youngTable$D4_12 <- NULL

for (i in c(1:length(youngClusters))) { assign(paste0("yng",youngClusters[i]), t(subset(youngTable, row.names(youngTable) %in% youngClusters[i]))) }


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




makeStateFig <- function (cptree, df, name="name")
{
    GSgradient <- contMap(cptree, df, res=1000, plot=FALSE, lwd=0.5, fsize=1, sig=0)
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
    GSgradient$tree$tip.label <- gsub("D10_7",paste0("G.turnerii 7, ",round(df[["D10_7"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D10_8",paste0("G.turnerii 8, ",round(df[["D10_8"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_2",paste0("G.harknessii, ",round(df[["D2_2"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D2_1_6",paste0("G.armourianum 6, ",round(df[["D2_1_6"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D11_1",paste0("G.schwendimanii 1, ",round(df[["D11_1"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_157",paste0("G.lobatum 157, ",round(df[["D7_157"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D7_4",paste0("G.lobatum 7, ",round(df[["D7_4"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D4_185",paste0("G.aridum 185, ",round(df[["D4_185"]])),GSgradient$tree$tip.label)
    GSgradient$tree$tip.label <- gsub("D9_4",paste0("G.laxum 4, ",round(df[["D9_4"]])),GSgradient$tree$tip.label)
    GSfit <- fastAnc(cptree, df, vars=TRUE, CI=TRUE)
    gfit <- round(GSfit$ace)
    
    png(paste0(name,".anc.png"), 5000, 3750, pointsize=12, res=600)
    plot(GSgradient, legend=FALSE, lwd=c(2,1), outline=FALSE, fsize=0.8) + nodelabels(round(gfit), adj=c(-0.3,0.2), frame="none", cex=0.6)
    dev.off()

    assign(paste0("G",name), GSgradient, envir=globalenv())
    assign(paste0("fit",name),gfit,envir=globalenv())
}

allTable <- Kbamount[,-c(1:9)]
names(allTable) <- gsub("\\.", "_", names(allTable))
allTable$D4_12 <- NULL




for (i in c(1:nrow(allTable))) { assign(paste0("all",row.names(allTable)[i]), t(allTable[i,])) }

allNames <- grep("allCL", ls(), value=TRUE)

for (cop in allNames) {  
	obj <- get(cop)
	names(obj) <- row.names(obj)
	assign(cop, obj, envir=globalenv())
}

cptree <- drop.tip(cptree, "D10_3")
for (name in allNames) { makeStateFig(cptree, get(name),(paste0("anc",name))) }



             
                             
                             

