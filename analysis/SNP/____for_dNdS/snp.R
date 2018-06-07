

# grep "#CHROM\|^Chr" all.PASS_snp.vcf > all.noHeader.snp.vcf
# for i in {5..18}; do awk '{ sub(/:.*$/, "", $i) } $0 ' all.oneAccession.vcf > temp.vcf; mv temp.vcf > all.oneAccession.vcf; done
# sed 's|\./\.|NA|g' all.oneAccession.vcf | sed 's|1/1|1|g' | sed 's|0/0|0|g' | sed 's|0/[123]|NA|g' | sed 's|2/2|2|g' | sed 's|3/3|3|g' | sed 's|1/[23]|NA|g' | sed 's|2/3|NA|g' > final.vcf
# grep -v NA final.vcf > final.noNA.vcf



#### Begin R
library(tidyr)
library(plyr)

## nt changes relative to F1-outgroup
snps <- read.table("final.noNA.vcf", header=T, stringsAsFactors=F)

	
Fdis <- data.frame(D1.35=double(),D10.7=double(),D11.1=double(),D2.1.6=double(),D2.2=double(),D3D.27=double(),D3K.57=double(),D4.185=double(),D5.8=double(),D6.5=double(),D7.157=double(),D8.8=double(),D9.4=double(),F1.1=double(), stringsAsFactors=FALSE)

for (D in names(snps)[c(5:18)]) {
    Fdis["distF",D] <- length(which(snps[,D] != snps[,18]))/nrow(snps)
}

write.table(Fdis, file="snp.distance", quote=F, sep="\t")

numSNP <- data.frame(D1.35=double(),D10.7=double(),D11.1=double(),D2.1.6=double(),D2.2=double(),D3D.27=double(),D3K.57=double(),D4.185=double(),D5.8=double(),D6.5=double(),D7.157=double(),D8.8=double(),D9.4=double(),F1.1=double(), stringsAsFactors=FALSE)

for (D in names(snps)[c(5:18)]) {
    numSNP["distF",D] <- length(which(snps[,D] != snps[,18]))
}

write.table(numSNP, file="snp.number", quote=F, sep="\t")




chrSNP <- data.frame(D1.35=double(),D10.7=double(),D11.1=double(),D2.1.6=double(),D2.2=double(),D3D.27=double(),D3K.57=double(),D4.185=double(),D5.8=double(),D6.5=double(),D7.157=double(),D8.8=double(),D9.4=double(),F1.1=double(), stringsAsFactors=FALSE)

for (D in names(snps)[c(5:18)]) {
    for (chr in (unique(snps$CHROM))) {
        chrSNP[chr,D] <- length(which(snps[snps$CHROM==chr,D] != snps[snps$CHROM==chr,18]))

}}

write.table(chrSNP, file="snp.number.chr", quote=F, sep="\t")