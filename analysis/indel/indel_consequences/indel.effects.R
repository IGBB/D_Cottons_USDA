indels <- read.table("D.snpeff.indel.filtered.annotations", header=T, sep="\t")
names(indels) <- c("CHROM","POS","REF","ALT","GENE","D1","D10","D11","D2_1","D2_2","D3D","D3K","D4","D5","D6","D7","D8","D9","EFF0","EFF1","EFF2","EFF3","EFF4","EFF5")
effTable <- t(data.frame(as.list(table(indels$EFF0))))

##### look at effects, broadly
100*sum(effTable[grep("^conservative", row.names(effTable))])/nrow(indels)
# [1] 18.97884

100*sum(effTable[grep("^frameshift_variant$", row.names(effTable))])/nrow(indels)
# [1] 48.85975

100*(1-sum(effTable[grep("^frameshift_variant$", row.names(effTable))])/sum(effTable[grep("^frameshift*", row.names(effTable))]))
# [1] 7.827438

100*sum(effTable[grep("disruptive", row.names(effTable))])/nrow(indels)
# [1] 27.20765

sum(effTable[grep("disruptive_inframe_deletion", row.names(effTable))])/sum(effTable[grep("disruptive_inframe_insertion", row.names(effTable))])
# [1] 2.21482

100*sum(effTable[grep("[.]disruptive_inframe_", row.names(effTable))])/sum(effTable[grep("disruptive", row.names(effTable))])
# [1] 1.396973

100*sum(effTable[grep("fusion", row.names(effTable))])/nrow(indels)
# [1] 0.1013556


#####################

length(unique(as.character(indels$GENE)))
# [1] 9341


length(unique(as.character(indels[-(grep("conservative", indels$EFF0)),]$GENE)))
# [1] 7767
gtLen <- indels[-(grep("conservative", indels$EFF0)),]

species <- names(indels[,c(6:18)])

nonref <- c()
geneind <- c()
gtLenGene <- c()
for (indiv in species) { 
x <- nrow(indels[indels[,indiv]>0,])
y <- length(unique(as.character(indels[indels[,indiv]>0,]$GENE)))
z <- length(unique(as.character(gtLen [gtLen [,indiv]>0,]$GENE)))
nonref <- c(nonref,x)
geneind <- c(geneind,y)
gtLenGene <- c(gtLenGene ,z)
}

spindels <- data.frame(species=species,indels=nonref,genes=geneind, gtLen=gtLenGene)
spNoD5 <- spindels[-9,]

mean(spNoD5$genes)
# [1] 2716.667

mean(spNoD5$gtLen)
