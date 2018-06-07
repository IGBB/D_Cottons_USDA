#### grep '^chr6' <file> | cut -f 1,$(seq -s, 6 13 <number of columns>)
# tail -c +2 all.PASS_indel.noHeader.vcf | head -n1 | cut -f 1,2,4,5,$(seq -s, 10 13 900) > Indel.header.final
# grep '^Chr01' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 10 13 900) >> Indel.header.final
# grep '^Chr02' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 15 13 900) >> Indel.header.final
# grep '^Chr03' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 16 13 900) >> Indel.header.final
# grep '^Chr04' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 17 13 900) >> Indel.header.final
# grep '^Chr05' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 18 13 900) >> Indel.header.final
# grep '^Chr06' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 19 13 900) >> Indel.header.final
# grep '^Chr07' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 20 13 900) >> Indel.header.final
# grep '^Chr08' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 21 13 900) >> Indel.header.final
# grep '^Chr09' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 22 13 900) >> Indel.header.final
# grep '^Chr10' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 11 13 900) >> Indel.header.final
# grep '^Chr11' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 12 13 900) >> Indel.header.final
# grep '^Chr12' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 13 13 900) >> Indel.header.final
# grep '^Chr13' all.PASS_indel.noHeader.vcf | cut -f 1,2,4,5,$(seq -s, 14 13 900) >> Indel.header.final

#### Begin R
library(tidyr)
library(plyr)


# read in all lines, transform designations
df <- read.table("Indel.header.final", header=T)
df <- df[,-c(5:14,19,23)]
names(df) <- gsub("\\.variant", "", names(df))
names(df) <- gsub("D3\\.", "D3", names(df))
names(df) <- gsub("\\.G.*", "", names(df))

df <- data.frame(lapply(df, function(x) { gsub("[.][/][.]", "NA", x, perl=TRUE) } ))
df <- data.frame(lapply(df, function(x) { gsub("\\d[/]\\d:\\d$", "NA", x, perl=TRUE) } ))

df <- data.frame(lapply(df, function(x) { gsub("0[/]0:.*", "0", x, perl=TRUE) } ))
df <- data.frame(lapply(df, function(x) { gsub("1[/]1:.*", "1", x, perl=TRUE) } ))
df <- data.frame(lapply(df, function(x) { gsub("2[/]2:.*", "2", x, perl=TRUE) } ))
df <- data.frame(lapply(df, function(x) { gsub("3[/]3:.*", "3", x, perl=TRUE) } ))
df <- data.frame(lapply(df, function(x) { gsub("4[/]4:.*", "4", x, perl=TRUE) } ))
df <- data.frame(lapply(df, function(x) { gsub("5[/]5:.*", "5", x, perl=TRUE) } ))
df <- data.frame(lapply(df, function(x) { gsub("6[/]6:.*", "6", x, perl=TRUE) } ))

df <- data.frame(lapply(df, function(x) { gsub("\\d[/]\\d:.*", "NA", x, perl=TRUE) } ))

df[df == "NA"] <- NA

save(df, file="pruned.indel.table.halfNA.Rdata")

keep <- apply(df[,c(5:32)], 1, function(x) length(unique(x[!is.na(x)])) != 1)
df <- df[keep, ]

nrow(df)
# [1] 7,174,698

df.na <- data.frame(NAs=integer(), stringsAsFactors=FALSE)
for (name in names(df[5:32])) {
    df.na[name,1] <- table(is.na(df[,name]))[[2]]
}

Findels <- df[,c("CHROM", "POS", "REF", "ALT","D1.35","D10.7","D11.1","D2.1.6","D2.2","D3D.27","D3K.57","D4.185","D5.8","D6.5","D7.157","D8.8","D9.4", "F1.1")]

indels <- Findels
indels$F1.1 <- NULL

Findels$NAsum <- rowSums(is.na(Findels))
indels$NAsum <- rowSums(is.na(indels))

Findels <- Findels[Findels$NAsum < 1, ]
nrow(Findels)
# [1] 1149943

indels <- indels[indels$NAsum < 1, ]
nrow(indels)
# [1] 1459578

save(Findels, file="pruned.Findel.table.noNA.Rdata")
save(indels, file="pruned.indel.table.noNA.Rdata")

indels$NAsum <- NULL
Findels$NAsum <- NULL

#################################
### write out the 
for (num in c(5:18) ) {
    filename <- names(Findels[num])
    write(as.character(Findels[,num]), file=filename)
	}
	

### make a tree with raxml
	
# to get these coded indels into phylip format, use this command in bash
# for a in [DF]*; do name=`(ls $a | cut -c 1-10 | xargs printf "%-10s\n" | tr ' ' 'Q' )`; cat $a | tr -d "\n" | fold -w50 | paste -sd'Q' - | sed 's/Q/QQQQQQQQQQ/g' | cat <(echo $name) - | tr -d "\n" | fold -w10 | paste -sd " " - | fold -w66 | paste -sd "\n" - | sed 's/Q/ /g' > $a.phy; done && numfiles=`ls [DF]*.phy | wc -l` && numchar=`wc -l D10.3 | cut -f1 -d ' '` && firstline=`echo $numfiles $numchar` && paste -d '\n' /dev/null *.phy | tail -n +2 | cat <( echo $firstline) - > all.indel.phy
#################################

indels <- separate(indels, 'ALT', paste("ALT", 1:6, sep="."), sep =",", extra="drop")
indels[is.na(indels)] <-""

indels$refLen<- nchar(as.vector(indels$REF))
indels$altLen.1<- nchar(as.vector(indels$ALT.1))
indels$altLen.2<- nchar(as.vector(indels$ALT.2))
indels$altLen.3<- nchar(as.vector(indels$ALT.3))
indels$altLen.4<- nchar(as.vector(indels$ALT.4))
indels$altLen.5<- nchar(as.vector(indels$ALT.5))
indels$altLen.6<- nchar(as.vector(indels$ALT.6))


for (col in c(10:34)) {
    indels[,col] <- as.numeric(as.character(indels[,col]))
}

keep <- apply(indels[,c(10:34)], 1, function(x) length(unique(x[!is.na(x)])) != 1)
indels <- indels[keep, ]

save(indels, file="indel.df.Rdata")

for (col in c(10:34)) {
    indels[,col] <- as.numeric(as.character(indels[,col]))
}

for (a in c(1:nrow(indels))) {
    for (i in c(10:34) ) {
	    indels[a,i] <- { ifelse(indels[a,i]==0, indels$refLen[a], 
		   ifelse(indels[a,i]==1, indels$altLen.1[a], 
		 	ifelse(indels[a,i]==2, indels$altLen.2[a],
				ifelse(indels[a,i]==3, indels$altLen.3[a], 
					ifelse(indels[a,i]==4, indels$altLen.4[a], 
						ifelse(indels[a,i]==5, indels$altLen.5[a], 
						    ifelse(indels[a,i]==6, indels$altLen.6[a], "err")))))))
	}
}}

save(indels, file="indel.df.replace.Rdata")


indels.relref <- indels[,c(1:2,10:34)]

for (i in c(3:27) ) {
	indels.relref[,i] <- indels.relref[,i]-indels[,35]
}


indels.accessions <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(indels.relref[,3:27])) {
      indels.accessions[i,"Chr01"] <- nrow(subset(indels.relref, CHROM=="Chr01" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr02"] <- nrow(subset(indels.relref, CHROM=="Chr02" & abs(indels.relref[,i])>0)) 
	indels.accessions[i,"Chr03"] <- nrow(subset(indels.relref, CHROM=="Chr03" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr04"] <- nrow(subset(indels.relref, CHROM=="Chr04" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr05"] <- nrow(subset(indels.relref, CHROM=="Chr05" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr06"] <- nrow(subset(indels.relref, CHROM=="Chr06" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr07"] <- nrow(subset(indels.relref, CHROM=="Chr07" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr08"] <- nrow(subset(indels.relref, CHROM=="Chr08" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr09"] <- nrow(subset(indels.relref, CHROM=="Chr09" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr10"] <- nrow(subset(indels.relref, CHROM=="Chr10" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr11"] <- nrow(subset(indels.relref, CHROM=="Chr11" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr12"] <- nrow(subset(indels.relref, CHROM=="Chr12" & abs(indels.relref[,i])>0))
	indels.accessions[i,"Chr13"] <- nrow(subset(indels.relref, CHROM=="Chr13" & abs(indels.relref[,i])>0)) }

indels.accessions$allChr <- rowSums(indels.accessions)

write.table(indels.accessions, file="NEW.indels.by.accession.tbl", quote=F, sep="\t")



indels.accessions.ins <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(indels.relref[,3:27])) {
      indels.accessions.ins[i,"Chr01"] <- nrow(subset(indels.relref, CHROM=="Chr01" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr02"] <- nrow(subset(indels.relref, CHROM=="Chr02" & indels.relref[,i]>0)) 
	indels.accessions.ins[i,"Chr03"] <- nrow(subset(indels.relref, CHROM=="Chr03" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr04"] <- nrow(subset(indels.relref, CHROM=="Chr04" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr05"] <- nrow(subset(indels.relref, CHROM=="Chr05" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr06"] <- nrow(subset(indels.relref, CHROM=="Chr06" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr07"] <- nrow(subset(indels.relref, CHROM=="Chr07" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr08"] <- nrow(subset(indels.relref, CHROM=="Chr08" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr09"] <- nrow(subset(indels.relref, CHROM=="Chr09" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr10"] <- nrow(subset(indels.relref, CHROM=="Chr10" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr11"] <- nrow(subset(indels.relref, CHROM=="Chr11" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr12"] <- nrow(subset(indels.relref, CHROM=="Chr12" & indels.relref[,i]>0))
	indels.accessions.ins[i,"Chr13"] <- nrow(subset(indels.relref, CHROM=="Chr13" & indels.relref[,i]>0)) }

indels.accessions.ins$allChr <- rowSums(indels.accessions.ins)

write.table(indels.accessions.ins, file="NEW.indels.ins.by.accession.tbl", quote=F, sep="\t")



indels.accessions.del <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(indels.relref[,3:27])) {
      indels.accessions.del[i,"Chr01"] <- nrow(subset(indels.relref, CHROM=="Chr01" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr02"] <- nrow(subset(indels.relref, CHROM=="Chr02" & indels.relref[,i]<0)) 
	indels.accessions.del[i,"Chr03"] <- nrow(subset(indels.relref, CHROM=="Chr03" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr04"] <- nrow(subset(indels.relref, CHROM=="Chr04" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr05"] <- nrow(subset(indels.relref, CHROM=="Chr05" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr06"] <- nrow(subset(indels.relref, CHROM=="Chr06" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr07"] <- nrow(subset(indels.relref, CHROM=="Chr07" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr08"] <- nrow(subset(indels.relref, CHROM=="Chr08" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr09"] <- nrow(subset(indels.relref, CHROM=="Chr09" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr10"] <- nrow(subset(indels.relref, CHROM=="Chr10" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr11"] <- nrow(subset(indels.relref, CHROM=="Chr11" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr12"] <- nrow(subset(indels.relref, CHROM=="Chr12" & indels.relref[,i]<0))
	indels.accessions.del[i,"Chr13"] <- nrow(subset(indels.relref, CHROM=="Chr13" & indels.relref[,i]<0)) }

indels.accessions.del$allChr <- rowSums(indels.accessions.del)

write.table(indels.accessions.del, file="NEW.indels.del.by.accession.tbl", quote=F, sep="\t")

indels.chr <- data.frame(indels=integer(), stringsAsFactors=F)

for (name in as.character(unique(indels$CHROM))) {
     indels.chr[name,"indels"] <- nrow(indels.relref[which(indels.relref$CHROM==name),])
}

write.table(indels.chr, file="NEW.indels.positions.chr.tbl", quote=F, sep="\t")


indels.relref.sm <- indels.relref[,c("CHROM", "POS","D1.35","D10.7","D11.1","D2.1.6","D2.2","D3D.27","D3K.57","D4.185","D5.8","D6.5","D7.157","D8.8","D9.4")]
keep.sm <- apply(indels.relref.sm[,c(3:15)], 1, function(x) length(unique(x[!is.na(x)])) != 1)
indels.relref.sm <- indels.relref.sm[keep.sm, ]


indels.chr.sm <- data.frame(indels=integer(), stringsAsFactors=F)

for (name in as.character(unique(indels$CHROM))) {
     indels.chr.sm[name,"indels"] <- nrow(indels.relref.sm[which(indels.relref.sm$CHROM==name),])
}

write.table(indels.chr.sm, file="Small.indels.positions.chr.tbl", quote=F, sep="\t")

########## F indels

Findels <- separate(Findels, 'ALT', paste("ALT", 1:6, sep="."), sep =",", extra="drop")
Findels[is.na(Findels)] <-""

Findels$refLen<- nchar(as.vector(Findels$REF))
Findels$altLen.1<- nchar(as.vector(Findels$ALT.1))
Findels$altLen.2<- nchar(as.vector(Findels$ALT.2))
Findels$altLen.3<- nchar(as.vector(Findels$ALT.3))
Findels$altLen.4<- nchar(as.vector(Findels$ALT.4))
Findels$altLen.5<- nchar(as.vector(Findels$ALT.5))
Findels$altLen.6<- nchar(as.vector(Findels$ALT.6))


for (col in c(10:23)) {
    Findels[,col] <- as.numeric(as.character(Findels[,col]))
}

save(Findels, file="Findel.df.Rdata")

for (a in c(1:nrow(Findels))) {
    for (i in c(10:23) ) {
	    Findels[a,i] <- { ifelse(Findels[a,i]==0, Findels$refLen[a], 
		   ifelse(Findels[a,i]==1, Findels$altLen.1[a], 
		 	ifelse(Findels[a,i]==2, Findels$altLen.2[a],
				ifelse(Findels[a,i]==3, Findels$altLen.3[a], 
					ifelse(Findels[a,i]==4, Findels$altLen.4[a], 
						ifelse(Findels[a,i]==5, Findels$altLen.5[a], 
						    ifelse(Findels[a,i]==6, Findels$altLen.6[a], "err")))))))
	}
}}

save(Findels, file="Findel.df.replace.Rdata")


Findels.relref <- Findels[,c(1:2,10:23)]

for (i in c(3:16) ) {
	Findels.relref[,i] <- Findels.relref[,i]-Findels.relref[,16]
}


Findels.accessions <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Findels.relref[,3:16])) {
      Findels.accessions[i,"Chr01"] <- nrow(subset(Findels.relref, CHROM=="Chr01" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr02"] <- nrow(subset(Findels.relref, CHROM=="Chr02" & abs(Findels.relref[,i])>0)) 
	Findels.accessions[i,"Chr03"] <- nrow(subset(Findels.relref, CHROM=="Chr03" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr04"] <- nrow(subset(Findels.relref, CHROM=="Chr04" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr05"] <- nrow(subset(Findels.relref, CHROM=="Chr05" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr06"] <- nrow(subset(Findels.relref, CHROM=="Chr06" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr07"] <- nrow(subset(Findels.relref, CHROM=="Chr07" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr08"] <- nrow(subset(Findels.relref, CHROM=="Chr08" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr09"] <- nrow(subset(Findels.relref, CHROM=="Chr09" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr10"] <- nrow(subset(Findels.relref, CHROM=="Chr10" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr11"] <- nrow(subset(Findels.relref, CHROM=="Chr11" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr12"] <- nrow(subset(Findels.relref, CHROM=="Chr12" & abs(Findels.relref[,i])>0))
	Findels.accessions[i,"Chr13"] <- nrow(subset(Findels.relref, CHROM=="Chr13" & abs(Findels.relref[,i])>0)) 
}

Findels.accessions$allChr <- rowSums(Findels.accessions)

write.table(Findels.accessions, file="Findels.by.accession.tbl", quote=F, sep="\t")



Findels.accessions.ins <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Findels.relref[,3:16])) {
      Findels.accessions.ins[i,"Chr01"] <- nrow(subset(Findels.relref, CHROM=="Chr01" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr02"] <- nrow(subset(Findels.relref, CHROM=="Chr02" & Findels.relref[,i]>0)) 
	Findels.accessions.ins[i,"Chr03"] <- nrow(subset(Findels.relref, CHROM=="Chr03" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr04"] <- nrow(subset(Findels.relref, CHROM=="Chr04" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr05"] <- nrow(subset(Findels.relref, CHROM=="Chr05" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr06"] <- nrow(subset(Findels.relref, CHROM=="Chr06" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr07"] <- nrow(subset(Findels.relref, CHROM=="Chr07" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr08"] <- nrow(subset(Findels.relref, CHROM=="Chr08" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr09"] <- nrow(subset(Findels.relref, CHROM=="Chr09" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr10"] <- nrow(subset(Findels.relref, CHROM=="Chr10" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr11"] <- nrow(subset(Findels.relref, CHROM=="Chr11" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr12"] <- nrow(subset(Findels.relref, CHROM=="Chr12" & Findels.relref[,i]>0))
	Findels.accessions.ins[i,"Chr13"] <- nrow(subset(Findels.relref, CHROM=="Chr13" & Findels.relref[,i]>0)) }

Findels.accessions.ins$allChr <- rowSums(Findels.accessions.ins)

write.table(Findels.accessions.ins, file="F.indels.ins.by.accession.tbl", quote=F, sep="\t")



Findels.accessions.del <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Findels.relref[,3:16])) {
      Findels.accessions.del[i,"Chr01"] <- nrow(subset(Findels.relref, CHROM=="Chr01" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr02"] <- nrow(subset(Findels.relref, CHROM=="Chr02" & Findels.relref[,i]<0)) 
	Findels.accessions.del[i,"Chr03"] <- nrow(subset(Findels.relref, CHROM=="Chr03" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr04"] <- nrow(subset(Findels.relref, CHROM=="Chr04" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr05"] <- nrow(subset(Findels.relref, CHROM=="Chr05" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr06"] <- nrow(subset(Findels.relref, CHROM=="Chr06" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr07"] <- nrow(subset(Findels.relref, CHROM=="Chr07" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr08"] <- nrow(subset(Findels.relref, CHROM=="Chr08" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr09"] <- nrow(subset(Findels.relref, CHROM=="Chr09" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr10"] <- nrow(subset(Findels.relref, CHROM=="Chr10" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr11"] <- nrow(subset(Findels.relref, CHROM=="Chr11" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr12"] <- nrow(subset(Findels.relref, CHROM=="Chr12" & Findels.relref[,i]<0))
	Findels.accessions.del[i,"Chr13"] <- nrow(subset(Findels.relref, CHROM=="Chr13" & Findels.relref[,i]<0)) }

Findels.accessions.del$allChr <- rowSums(Findels.accessions.del)

write.table(Findels.accessions.del, file="Findels.del.by.accession.tbl", quote=F, sep="\t")

Findels.minmax <- Findels[,24:30]

Findels.minmax[(Findels.minmax==0)] <-NA

Findels.minmax$size <- apply(Findels.minmax,1,max, na.rm=T) - apply(Findels.minmax,1,min, na.rm=T)
Findels.minmax$min <- apply(Findels.minmax,1,min, na.rm=T)
Findels.minmax$max <- apply(Findels.minmax,1,max, na.rm=T) 

Findels.minmax$CHROM <- Findels$CHROM

Findels.chr <- data.frame(indels=integer(), avgSize=integer(), min=integer(), max=integer(), stringsAsFactors=F)

for (name in as.character(unique(Findels.minmax$CHROM))) {
     Findels.chr[name,"indels"] <- nrow(Findels.minmax[which(Findels.minmax$CHROM==name),])
     Findels.chr[name,"avgSize"] <- mean(Findels.minmax[which(Findels.minmax$CHROM==name),]$size)
     Findels.chr[name,"min"] <- min(Findels.minmax[which(Findels.minmax$CHROM==name),]$min)
     Findels.chr[name,"max"] <- max(Findels.minmax[which(Findels.minmax$CHROM==name),]$max)

}

write.table(Findels.chr, file="Findels.positions.chr.tbl", quote=F, sep="\t")

save(Findels.chr, Findels.accessions.ins, Findels.accessions.del, Findels.accessions, Findels, Findels.relref, file="Findels.Rdata")

####### now Findels that occured within the Ds ######


Dkeep <- apply(Findels[,c(10:22)], 1, function(x) length(unique(x[!is.na(x)])) != 1)
Dindels <- Findels[Dkeep, ]

nrow(Dindels)
# [1] 761746


Dindels.relref <- Dindels[,c(1:2,10:23)]

for (i in c(3:16) ) {
	Dindels.relref[,i] <- Dindels.relref[,i]-Dindels.relref[,16]
}


Dindels.accessions <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Dindels.relref[,3:16])) {
      Dindels.accessions[i,"Chr01"] <- nrow(subset(Dindels.relref, CHROM=="Chr01" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr02"] <- nrow(subset(Dindels.relref, CHROM=="Chr02" & abs(Dindels.relref[,i])>0)) 
	Dindels.accessions[i,"Chr03"] <- nrow(subset(Dindels.relref, CHROM=="Chr03" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr04"] <- nrow(subset(Dindels.relref, CHROM=="Chr04" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr05"] <- nrow(subset(Dindels.relref, CHROM=="Chr05" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr06"] <- nrow(subset(Dindels.relref, CHROM=="Chr06" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr07"] <- nrow(subset(Dindels.relref, CHROM=="Chr07" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr08"] <- nrow(subset(Dindels.relref, CHROM=="Chr08" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr09"] <- nrow(subset(Dindels.relref, CHROM=="Chr09" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr10"] <- nrow(subset(Dindels.relref, CHROM=="Chr10" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr11"] <- nrow(subset(Dindels.relref, CHROM=="Chr11" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr12"] <- nrow(subset(Dindels.relref, CHROM=="Chr12" & abs(Dindels.relref[,i])>0))
	Dindels.accessions[i,"Chr13"] <- nrow(subset(Dindels.relref, CHROM=="Chr13" & abs(Dindels.relref[,i])>0)) 
}

Dindels.accessions$allChr <- rowSums(Dindels.accessions)

write.table(Dindels.accessions, file="Dindels.by.accession.tbl", quote=F, sep="\t")



Dindels.accessions.ins <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Dindels.relref[,3:16])) {
      Dindels.accessions.ins[i,"Chr01"] <- nrow(subset(Dindels.relref, CHROM=="Chr01" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr02"] <- nrow(subset(Dindels.relref, CHROM=="Chr02" & Dindels.relref[,i]>0)) 
	Dindels.accessions.ins[i,"Chr03"] <- nrow(subset(Dindels.relref, CHROM=="Chr03" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr04"] <- nrow(subset(Dindels.relref, CHROM=="Chr04" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr05"] <- nrow(subset(Dindels.relref, CHROM=="Chr05" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr06"] <- nrow(subset(Dindels.relref, CHROM=="Chr06" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr07"] <- nrow(subset(Dindels.relref, CHROM=="Chr07" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr08"] <- nrow(subset(Dindels.relref, CHROM=="Chr08" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr09"] <- nrow(subset(Dindels.relref, CHROM=="Chr09" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr10"] <- nrow(subset(Dindels.relref, CHROM=="Chr10" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr11"] <- nrow(subset(Dindels.relref, CHROM=="Chr11" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr12"] <- nrow(subset(Dindels.relref, CHROM=="Chr12" & Dindels.relref[,i]>0))
	Dindels.accessions.ins[i,"Chr13"] <- nrow(subset(Dindels.relref, CHROM=="Chr13" & Dindels.relref[,i]>0)) }

Dindels.accessions.ins$allChr <- rowSums(Dindels.accessions.ins)

write.table(Dindels.accessions.ins, file="Dindels.ins.by.accession.tbl", quote=F, sep="\t")



Dindels.accessions.del <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Dindels.relref[,3:16])) {
      Dindels.accessions.del[i,"Chr01"] <- nrow(subset(Dindels.relref, CHROM=="Chr01" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr02"] <- nrow(subset(Dindels.relref, CHROM=="Chr02" & Dindels.relref[,i]<0)) 
	Dindels.accessions.del[i,"Chr03"] <- nrow(subset(Dindels.relref, CHROM=="Chr03" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr04"] <- nrow(subset(Dindels.relref, CHROM=="Chr04" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr05"] <- nrow(subset(Dindels.relref, CHROM=="Chr05" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr06"] <- nrow(subset(Dindels.relref, CHROM=="Chr06" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr07"] <- nrow(subset(Dindels.relref, CHROM=="Chr07" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr08"] <- nrow(subset(Dindels.relref, CHROM=="Chr08" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr09"] <- nrow(subset(Dindels.relref, CHROM=="Chr09" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr10"] <- nrow(subset(Dindels.relref, CHROM=="Chr10" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr11"] <- nrow(subset(Dindels.relref, CHROM=="Chr11" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr12"] <- nrow(subset(Dindels.relref, CHROM=="Chr12" & Dindels.relref[,i]<0))
	Dindels.accessions.del[i,"Chr13"] <- nrow(subset(Dindels.relref, CHROM=="Chr13" & Dindels.relref[,i]<0)) }

Dindels.accessions.del$allChr <- rowSums(Dindels.accessions.del)

write.table(Dindels.accessions.del, file="Dindels.del.by.accession.tbl", quote=F, sep="\t")

Dindels.minmax <- Dindels[,24:30]

Dindels.minmax[(Dindels.minmax==0)] <-NA

Dindels.minmax$size <- apply(Dindels.minmax,1,max, na.rm=T) - apply(Dindels.minmax,1,min, na.rm=T)
Dindels.minmax$min <- apply(Dindels.minmax,1,min, na.rm=T)
Dindels.minmax$max <- apply(Dindels.minmax,1,max, na.rm=T) 

Dindels.minmax$CHROM <- Dindels$CHROM

Dindels.chr <- data.frame(indels=integer(), avgSize=integer(), min=integer(), max=integer(), stringsAsFactors=F)

for (name in as.character(unique(Dindels.minmax$CHROM))) {
     Dindels.chr[name,"indels"] <- nrow(Dindels.minmax[which(Dindels.minmax$CHROM==name),])
     Dindels.chr[name,"avgSize"] <- mean(Dindels.minmax[which(Dindels.minmax$CHROM==name),]$size)
     Dindels.chr[name,"min"] <- min(Dindels.minmax[which(Dindels.minmax$CHROM==name),]$min)
     Dindels.chr[name,"max"] <- max(Dindels.minmax[which(Dindels.minmax$CHROM==name),]$max)

}

write.table(Dindels.chr, file="Dindels.positions.chr.tbl", quote=F, sep="\t")

save(Dindels.chr, Dindels.accessions.ins, Dindels.accessions.del, Dindels.accessions, Dindels, Dindels.relref, file="Dindels.Rdata")


















Dindels.accessions.ins.len <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Dindels.relref[,3:16])) {
      Dindels.accessions.ins.len[i,"Chr01"] <- sum(subset(Dindels.relref, CHROM=="Chr01" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr02"] <- sum(subset(Dindels.relref, CHROM=="Chr02" & Dindels.relref[,i]>0)[,i]) 
	Dindels.accessions.ins.len[i,"Chr03"] <- sum(subset(Dindels.relref, CHROM=="Chr03" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr04"] <- sum(subset(Dindels.relref, CHROM=="Chr04" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr05"] <- sum(subset(Dindels.relref, CHROM=="Chr05" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr06"] <- sum(subset(Dindels.relref, CHROM=="Chr06" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr07"] <- sum(subset(Dindels.relref, CHROM=="Chr07" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr08"] <- sum(subset(Dindels.relref, CHROM=="Chr08" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr09"] <- sum(subset(Dindels.relref, CHROM=="Chr09" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr10"] <- sum(subset(Dindels.relref, CHROM=="Chr10" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr11"] <- sum(subset(Dindels.relref, CHROM=="Chr11" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr12"] <- sum(subset(Dindels.relref, CHROM=="Chr12" & Dindels.relref[,i]>0)[,i])
	Dindels.accessions.ins.len[i,"Chr13"] <- sum(subset(Dindels.relref, CHROM=="Chr13" & Dindels.relref[,i]>0)[,i]) }

Dindels.accessions.ins.len$allChr <- rowSums(Dindels.accessions.ins.len)

write.table(Dindels.accessions.ins.len, file="Dindels.ins.len.by.accession.tbl", quote=F, sep="\t")



Dindels.accessions.del.len <- data.frame(Chr01=integer(), Chr02=integer(), Chr03=integer(), Chr04=integer(), Chr05=integer(), Chr06=integer(), Chr07=integer(), Chr08=integer(), Chr09=integer(), Chr10=integer(), Chr11=integer(), Chr12=integer(), Chr13=integer(), stringsAsFactors=FALSE)

for (i in names(Dindels.relref[,3:16])) {
      Dindels.accessions.del.len[i,"Chr01"] <- sum(subset(Dindels.relref, CHROM=="Chr01" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr02"] <- sum(subset(Dindels.relref, CHROM=="Chr02" & Dindels.relref[,i]<0)[,i]) 
	Dindels.accessions.del.len[i,"Chr03"] <- sum(subset(Dindels.relref, CHROM=="Chr03" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr04"] <- sum(subset(Dindels.relref, CHROM=="Chr04" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr05"] <- sum(subset(Dindels.relref, CHROM=="Chr05" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr06"] <- sum(subset(Dindels.relref, CHROM=="Chr06" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr07"] <- sum(subset(Dindels.relref, CHROM=="Chr07" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr08"] <- sum(subset(Dindels.relref, CHROM=="Chr08" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr09"] <- sum(subset(Dindels.relref, CHROM=="Chr09" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr10"] <- sum(subset(Dindels.relref, CHROM=="Chr10" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr11"] <- sum(subset(Dindels.relref, CHROM=="Chr11" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr12"] <- sum(subset(Dindels.relref, CHROM=="Chr12" & Dindels.relref[,i]<0)[,i])
	Dindels.accessions.del.len[i,"Chr13"] <- sum(subset(Dindels.relref, CHROM=="Chr13" & Dindels.relref[,i]<0)[,i]) }

Dindels.accessions.del.len$allChr <- rowSums(Dindels.accessions.del.len)

write.table(Dindels.accessions.del.len, file="Dindels.del.len.by.accession.tbl", quote=F, sep="\t")
