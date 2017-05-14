library(ggplot2)
library(scales)
library(ggrepel)
library(factoextra)
library(reshape2)
library(gridExtra)
library(geomorph)
library(testit)

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

ord_table <- annot_clust[,(3:35)]
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
lapply(Mb_table[24:27], function(x) x/910),
lapply(Mb_table[28], function(x) x/929),
lapply(Mb_table[5], function(x) x/910),
lapply(Mb_table[6:7], function(x) x/880),
lapply(Mb_table[8:9], function(x) x/919),
lapply(Mb_table[10:16], function(x) x/880),
lapply(Mb_table[17:18], function(x) x/841), 
lapply(Mb_table[19:20], function(x) x/934),
lapply(Mb_table[21:22], function(x) x/851),
lapply(Mb_table[23], function(x) x/934),
lapply(Mb_table[29:30], function(x) x/1667),
lapply(Mb_table[31:32], function(x) x/1698),
lapply(Mb_table[33], function(x) x/1311) )

perc_table <- perc_table[,c(1:28)]
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
png("cotton.ordination.png", 5000, 5000, pointsize=12, res=600)
ggplot(cmdpoints, aes(x=V1, y=V2, color=subsection)) + geom_point(size=2) + xlab("PCoA component 1") + ylab("PCoA component 2") + geom_text_repel(aes(label=species))
dev.off()

### PCoA of log-transformed data ###
# evaluates how close the overall repetitive profiles are to one another

ord_table <- ord_table[,c(1:28)]
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
ggplot(logcmdpoints, aes(x=V1, y=V2, color=subsection)) + geom_point(size=2) + xlab("PCoA component 1") + ylab("PCoA component 2") + geom_text_repel(aes(label=species))
dev.off()

### took out Procrustea ANOVA because nothing appears different

########### characterize composition ###########

annot_clust$cluster <- NULL

# 9.5 multiplier represents # reads (x) * 95nt/read * 1 kb/1000nt * 100% = # reads * 9.5 = # Kb in entire genome for that class 
Kbamount <- data.frame(annot_clust[1], apply(annot_clust[2:34], 2, function (x) x*9.5))
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
KBsum$D10sum <- rowMeans(KBsum[,25:28])
KBsum$D11sum <-KBsum$D11_1
KBsum$A1sum  <- rowMeans(KBsum[,30:31])
KBsum$A2sum  <- rowMeans(KBsum[,32:33])
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
KBsum$D10min <- apply(KBsum[,25:28], 1, min)
KBsum$D11min <-KBsum$D11_1
KBsum$A1min  <- apply(KBsum[,30:31], 1, min)
KBsum$A2min  <- apply(KBsum[,32:33], 1, min)
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
KBsum$D10max <- apply(KBsum[,25:28], 1, max)
KBsum$D11max <-KBsum$D11_1
KBsum$A1max  <- apply(KBsum[,30:31], 1, max)
KBsum$A2max  <- apply(KBsum[,32:33], 1, max)
KBsum$F1max <-KBsum$F1_1

KBsum <- KBsum[,-(2:34)]
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

sum(KBsum$kirkii)/1000 # convert to Mb
# [1] 110.3615
sum(KBsum$kokia_)/1000 # convert to Mb
# [1] 109.4685

sum(KBsum$D1sum)/1000
sum(KBsum$D2.1sum)/1000
sum(KBsum$D2.2sum)/1000
sum(KBsum$D3Dsum)/1000
sum(KBsum$D3Ksum)/1000
sum(KBsum$D4sum)/1000
sum(KBsum$D5sum)/1000
sum(KBsum$D6sum)/1000
sum(KBsum$D7sum)/1000
sum(KBsum$D8sum)/1000
sum(KBsum$D9sum)/1000
sum(KBsum$D10sum)/1000
sum(KBsum$D11sum)/1000
sum(KBsum$A1sum)/1000
sum(KBsum$A2sum)/1000
sum(KBsum$F1sum)/1000




### chi2 test to determine clusters that are significantly different ###
##
##
##
# ok this is where we should get creative.  a billion t-tests won't tell us anything
# very interesting, but what will tell us something interesting is when something changes
# so, when the phylogeny is done, we can do sister species tests for differential cluster
# abundace. If the clusters are not differentially abundant, then merge the species for 
# that cluster. If the species are differentially abundant, then do ancestral state 
# reconstruction. Now use ancestral state as counts for that cluster?









  
chi_table <- annot_clust[,c(1,15:16)]
chi_table <- chi_table[!(rowSums(chi_table[,c(2:3)])==0),]

chi_table$p.value <- apply(chi_table[,c(2:3)], 1, function(x) chisq.test(x, simulate.p.value=TRUE, B=50000)$p.value)
chi_table$statistic <- apply(chi_table[,c(2:3)], 1, function(x) chisq.test(x)$statistic)
# note if you DO NOT use the p-value simulation
# there are 35 clusters with the warning "In chisq.test(x) : Chi-squared approximation may be incorrect"
# The "simulate.p.value=T" option (default value is FALSE) does the Monte Carlo simulation; from https://ww2.coastal.edu/kingw/statistics/R-tutorials/independ.html

chi_table$p.adjust <- p.adjust(chi_table$p.value, method="BH")

# > row.names(chi_table[(chi_table$p.adjust<0.001),])
# [1] "CL0002" "CL0005" "CL0050" "CL0084" "CL0096" "CL0097" "CL0101" "CL0105" "CL0107" "CL0110" "CL0116" "CL0117" "CL0119"
# [14] "CL0121" "CL0126" "CL0129" "CL0141" "CL0162" "CL0164" "CL0175" "CL0177" "CL0187" "CL0188"
# > row.names(chi_table[(chi_table$p.adjust<0.005),])
# [1] "CL0002" "CL0005" "CL0017" "CL0050" "CL0066" "CL0082" "CL0084" "CL0085" "CL0096" "CL0097" "CL0101" "CL0105" "CL0107"
# [14] "CL0110" "CL0116" "CL0117" "CL0119" "CL0121" "CL0126" "CL0129" "CL0136" "CL0141" "CL0158" "CL0162" "CL0164" "CL0175"
# [27] "CL0177" "CL0187" "CL0188" "CL0190" "CL0191" "CL0238" "CL0253" "CL0274"
# > row.names(chi_table[(chi_table$p.adjust<0.01),])
# [1] "CL0002" "CL0005" "CL0017" "CL0050" "CL0066" "CL0069" "CL0082" "CL0084" "CL0085" "CL0096" "CL0097" "CL0098" "CL0100"
# [14] "CL0101" "CL0105" "CL0107" "CL0110" "CL0116" "CL0117" "CL0118" "CL0119" "CL0121" "CL0126" "CL0129" "CL0136" "CL0141"
# [27] "CL0149" "CL0158" "CL0162" "CL0164" "CL0175" "CL0177" "CL0187" "CL0188" "CL0190" "CL0191" "CL0203" "CL0238" "CL0242"
# [40] "CL0253" "CL0274"p05 <- row.names(chi_table[(chi_table$p.adjust<0.05),])

p05 <- row.names(chi_table[(chi_table$p.adjust<0.05),])
# [1] "CL0002" "CL0005" "CL0010" "CL0017" "CL0050" "CL0060" "CL0066" "CL0069" "CL0082" "CL0084" "CL0085" "CL0096" "CL0097"
# [14] "CL0098" "CL0100" "CL0101" "CL0105" "CL0107" "CL0110" "CL0116" "CL0117" "CL0118" "CL0119" "CL0121" "CL0126" "CL0128"
# [27] "CL0129" "CL0131" "CL0136" "CL0141" "CL0147" "CL0149" "CL0150" "CL0153" "CL0158" "CL0161" "CL0162" "CL0164" "CL0168"
# [40] "CL0170" "CL0175" "CL0177" "CL0187" "CL0188" "CL0190" "CL0191" "CL0202" "CL0203" "CL0219" "CL0238" "CL0242" "CL0243"
# [53] "CL0253" "CL0271" "CL0274"

#  clusters with warnings if no p-value simulation
warnME <- NULL
for (i in c(1:188)) { assign("last.warning", NULL, envir = baseenv()); ifelse(has_warning(chisq.test(chi_table[i,c(2:3)])), warnME <- c(warnME,(row.names(chi_table)[i])),"") }
#  [1] "CL0001" "CL0004" "CL0013" "CL0018" "CL0019" "CL0024" "CL0031" "CL0034" "CL0037" "CL0044" "CL0054" "CL0055" "CL0064"
# [14] "CL0067" "CL0080" "CL0092" "CL0132" "CL0134" "CL0155" "CL0179" "CL0181" "CL0195" "CL0209" "CL0214" "CL0215" "CL0223"
# [27] "CL0233" "CL0234" "CL0243" "CL0246" "CL0252" "CL0259" "CL0264" "CL0267" "CL0269"

intersect(warnME, p05)
#  only CL243 overlaps between the warnings and sig diff clusters at p<0.05
                             
# make a table of only significant clusters and then learn about it 
sigtable <- chi_table[p05,c(1:3)]
sigtable$greater <- ifelse(sigtable$kokia_ > sigtable$kirkii, "kokia", ifelse(sigtable$kirkii > sigtable$kokia_, "kirkii", "same"))
table(sigtable$greater)
#kirkii  kokia 
#    21     34 

table(sigtable$Lineage)
#  *    LTR     LTR/Copia      LTR/Gypsy       non-LTR_retroposon        Unspecified 
#  2     11         6              34                  1                      1 

sum(sigtable$kirkii)*0.0095 # 0.0095 is the same factor used above, with the conversion to Mb
# [1] 70.4235
sum(sigtable$kokia_)*0.0095
# [1] 68.894
sum(sigtable$kirkii)*0.0095-sum(sigtable$kokia_)*0.0095
# [1] 1.5295
                            
                             
                             
########### relative aging of transposable elements ###########

# see TE_dating_histogram.pl for table generation

ages <- read.table("cluster.ages", header=T, row.names=1, sep="\t")
KdGkages <- read.table("KokKirk.cluster.ages", header=T, row.names=1, sep="\t")
KdGkages$age <- KdGkages$Category

#Category	age
#1	young
#3	old
#4	young
#5	old
#4*	old

KdGkages$age <- sub("3|5|4\\*","old",KdGkages$age)
KdGkages$age <- sub("1|4","young",KdGkages$age)

table(ages$age)
#  old young 
#  202    72 

table(KdGkages$age)

# old young 
# 127    49 


sigages <- merge(sigtable,ages,by="row.names")
table(sigages$age)
#  old young 
#   30    25 

sigGkKd <- merge(sigtable, KdGkages, by="row.names")
table(sigGkKd$age)

# old young 
#  31    24 

changed <- NULL
for (i in c(1:55)) { ifelse(identical(as.character(sigages$age[i]), as.character(sigGkKd$age[i])), "", changed <- c(changed,sigages$Row.names[i])) }
changed
# CL0097 young -> old

