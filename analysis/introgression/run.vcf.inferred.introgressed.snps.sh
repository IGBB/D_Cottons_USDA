########################################
# Assessing introgression in D4-Colima
#
# Past analyses indicate that there is likely D3 introgression in D4-Colima,
# due to the presence of D3 cytoplasm in D4-Colima.
# This relationship was also demonstrated via AFLPs in Alvarez 2006
# Whole chloroplast phylogenies constructed here recapitulate Integrifolia-like chloroplast in D4-Colima
# ABBA-BABA tests are not perfectly applicable here because D6 itself is multiply introgressant
# from both an African-like source and from G. raimondii

###
# to determine possible shared, derived SNPs
# prune vcf to only include both D4 accessions, D3D, and D6
# require D3D != D6 && D4C != D4
# determine how many SNPs D3D == D4C
# and how many SNPs D3D == D4

module load bedtools2/2.24.0
# note: bedtools2/2.26.0 appears to have limitations that cause the gff to fail

grep '^#CHROM\|^Chr' all.PASS_snp.vcf | awk ' BEGIN{OFS="\t";} { print $1,$2,$4,$5,$30,$33,$34,$42  } ' | sed 's|[.]/[.]:|Q/Q:|g' |  awk -v OFS="\t"  ' { for ( i= 5; i<=23; i++) { gsub(":.*", "", $i) } print }  ' | awk ' {  if ( ($5 != "Q/Q" && $6 != "Q/Q" && $7 != "Q/Q" && $8 != "Q/Q" && $6 != $7 && $5 != $8) || ($1 == "#CHROM")) { print } } ' | sed '/0\/1/d' | sed '/0\/2/d' | sed '/0\/3/d' | sed '/1\/2/d' | sed '/1\/3/d' | sed '/2\/3/d' > D4.D3D.D6-5.vcf

awk ' {  if ( ($5 == $6 && $7 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3D.D6-5.vcf > D4CeqD3D.vcf
awk ' {  if ( ($5 == $7 && $6 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3D.D6-5.vcf > D4eqD3D.vcf
for i in {01..13}; do echo Chr$i >> D4CeqD3D.counts; grep -c Chr$i D4CeqD3D.vcf >> D4CeqD3D.counts; done
for i in {01..13}; do echo Chr$i >> D4eqD3D.counts; grep -c Chr$i D4eqD3D.vcf >> D4eqD3D.counts; done

awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4CeqD3D.vcf | sed '/#/d' > D4CeqD3D.bed
awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4eqD3D.vcf | sed '/#/d' > D4eqD3D.bed

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3D.bed | wc -l 
#7843

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3D.bed | wc -l 
#7419

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3D.bed | cut -f9| sort | uniq | wc -l
#4808

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3D.bed | cut -f9| sort | uniq | wc -l
#4721

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3D.bed | cut -f9| sort | uniq > D4CeqD3D.genes
bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3D.bed | cut -f9| sort | uniq > D4eqD3D.genes

######## R code to look at the overlap in genes
## D4c <- read.table("D4CeqD3D.genes", fill=TRUE, stringsAsFactors=F)
## D4 <- read.table("D4eqD3D.genes", fill=TRUE, stringsAsFactors=F)
## length(intersect(D4$V1, D4c$V1))
## [1] 1542
########





###
# do we see the same when we use D3K
# prune vcf to only include both D4 accessions, D3K, and D6
# require D3K != D6 && D4C != D4
# determine how many SNPs D3K == D4C
# and how many SNPs D3K == D4

module load bedtools2/2.24.0
# note: bedtools2/2.26.0 appears to have limitations that cause the gff to fail

grep '^#CHROM\|^Chr' all.PASS_snp.vcf | awk ' BEGIN{OFS="\t";} { print $1,$2,$4,$5,$31,$33,$34,$42  } ' | sed 's|[.]/[.]:|Q/Q:|g' |  awk -v OFS="\t"  ' { for ( i= 5; i<=23; i++) { gsub(":.*", "", $i) } print }  ' | awk ' {  if ( ($5 != "Q/Q" && $6 != "Q/Q" && $7 != "Q/Q" && $8 != "Q/Q" && $6 != $7 && $5 != $8) || ($1 == "#CHROM")) { print } } ' | sed '/0\/1/d' | sed '/0\/2/d' | sed '/0\/3/d' | sed '/1\/2/d' | sed '/1\/3/d' | sed '/2\/3/d' > D4.D3K.D6-5.vcf

awk ' {  if ( ($5 == $6 && $7 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3K.D6-5.vcf > D4CeqD3K.vcf
awk ' {  if ( ($5 == $7 && $6 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3K.D6-5.vcf > D4eqD3K.vcf
for i in {01..13}; do echo Chr$i >> D4CeqD3K.counts; grep -c Chr$i D4CeqD3K.vcf >> D4CeqD3K.counts; done
for i in {01..13}; do echo Chr$i >> D4eqD3K.counts; grep -c Chr$i D4eqD3K.vcf >> D4eqD3K.counts; done

awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4CeqD3K.vcf | sed '/#/d' > D4CeqD3K.bed
awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4eqD3K.vcf | sed '/#/d' > D4eqD3K.bed

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3K.bed | wc -l
#7825

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3K.bed | wc -l
#7381

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3K.bed | cut -f9| sort | uniq | wc -l
#4804

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3K.bed | cut -f9| sort | uniq | wc -l
#4705

###
# to determine possible shared, derived SNPs
# prune vcf to only include both D4 accessions, D3D, and D6
# require D3D != F1 && D4C != D4
# determine how many SNPs F == D4C
# and how many SNPs F == D4

module load bedtools2/2.24.0
# note: bedtools2/2.26.0 appears to have limitations that cause the gff to fail

grep '^#CHROM\|^Chr' all.PASS_snp.vcf | awk ' BEGIN{OFS="\t";} { print $1,$2,$4,$5,$30,$33,$34,$49  } ' | sed 's|[.]/[.]:|Q/Q:|g' |  awk -v OFS="\t"  ' { for ( i= 5; i<=23; i++) { gsub(":.*", "", $i) } print }  ' | awk ' {  if ( ($5 != "Q/Q" && $6 != "Q/Q" && $7 != "Q/Q" && $8 != "Q/Q" && $6 != $7 && $5 != $8) || ($1 == "#CHROM")) { print } } ' | sed '/0\/1/d' | sed '/0\/2/d' | sed '/0\/3/d' | sed '/1\/2/d' | sed '/1\/3/d' | sed '/2\/3/d' > D4.D3D.F.vcf

awk ' {  if ( ($5 == $6 && $7 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3D.F.vcf > D4CeqD3D.F.vcf
awk ' {  if ( ($5 == $7 && $6 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3D.F.vcf > D4eqD3D.F.vcf
for i in {01..13}; do echo Chr$i >> D4CeqD3D.F.counts; grep -c Chr$i D4CeqD3D.F.vcf >> D4CeqD3D.F.counts; done
for i in {01..13}; do echo Chr$i >> D4eqD3D.F.counts; grep -c Chr$i D4eqD3D.F.vcf >> D4eqD3D.F.counts; done

awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4CeqD3D.F.vcf | sed '/#/d' > D4CeqD3D.F.bed
awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4eqD3D.F.vcf | sed '/#/d' > D4eqD3D.F.bed

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3D.F.bed | wc -l
#9883

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3D.F.bed | wc -l
#8469

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3D.F.bed | cut -f9| sort | uniq | wc -l
#6355

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3D.F.bed | cut -f9| sort | uniq | wc -l
#5522

###
# do we see the same when we use D3K
# prune vcf to only include both D4 accessions, D3K, and F
# require D3K != F && D4C != D4
# determine how many SNPs D3K == D4C
# and how many SNPs D3K == D4

module load bedtools2/2.24.0
# note: bedtools2/2.26.0 appears to have limitations that cause the gff to fail

grep '^#CHROM\|^Chr' all.PASS_snp.vcf | awk ' BEGIN{OFS="\t";} { print $1,$2,$4,$5,$31,$33,$34,$49  } ' | sed 's|[.]/[.]:|Q/Q:|g' |  awk -v OFS="\t"  ' { for ( i= 5; i<=23; i++) { gsub(":.*", "", $i) } print }  ' | awk ' {  if ( ($5 != "Q/Q" && $6 != "Q/Q" && $7 != "Q/Q" && $8 != "Q/Q" && $6 != $7 && $5 != $8) || ($1 == "#CHROM")) { print } } ' | sed '/0\/1/d' | sed '/0\/2/d' | sed '/0\/3/d' | sed '/1\/2/d' | sed '/1\/3/d' | sed '/2\/3/d' > D4.D3K.F.vcf

awk ' {  if ( ($5 == $6 && $7 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3K.F.vcf > D4CeqD3K.F.vcf
awk ' {  if ( ($5 == $7 && $6 == $8) || ($1 == "#CHROM") ) { print } } ' D4.D3K.F.vcf > D4eqD3K.F.vcf
for i in {01..13}; do echo Chr$i >> D4CeqD3K.F.counts; grep -c Chr$i D4CeqD3K.F.vcf >> D4CeqD3K.F.counts; done
for i in {01..13}; do echo Chr$i >> D4eqD3K.F.counts; grep -c Chr$i D4eqD3K.F.vcf >> D4eqD3K.F.counts; done

awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4CeqD3K.F.vcf | sed '/#/d' > D4CeqD3K.F.bed
awk ' BEGIN{OFS="\t";} { print $1, $2, $2 }' D4eqD3K.F.vcf | sed '/#/d' > D4eqD3K.F.bed

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3K.F.bed | wc -l
#9871

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3K.F.bed | wc -l
#8448

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4CeqD3K.F.bed | cut -f9| sort | uniq | wc -l
#6359

bedtools intersect -nonamecheck -wo -a D5.CDS_range.gff -b D4eqD3K.F.bed | cut -f9| sort | uniq | wc -l
#5508

