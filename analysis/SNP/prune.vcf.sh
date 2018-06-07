########################################
# Assessing introgression in D4-Colima
#
# Past analyses indicate that there is likely D3 introgression in D4-Colima,
# due to the presence of D3 cytoplasm in D4-Colima.
# This relationship was also demonstrated via AFLPs in Alvarez 2006
# Whole chloroplast phylogenies constructed here WHATDOTHEYSAY
# ABBA-BABA tests are not perfectly applicable here because D6 itself is multiply introgressant
# from both an African-like source and from G. raimondii

grep '^#CHROM\|^Chr' all.PASS_snp.vcf | awk ' BEGIN{OFS="\t";} { print $1,$2,$4,$5,$30,$33,$34,$42  } ' | sed 's|[.]/[.]:|Q/Q:|g' |  awk -v OFS="\t"  ' { for ( i= 5; i<=23; i++) { gsub(":.*", "", $i) } print }  ' | awk ' {  if ( ($5 != "Q/Q" && $6 != "Q/Q" && $7 != "Q/Q" && $8 != "Q/Q" && $6 != $7 && $5 != $8) || ($1 == "#CHROM")) { print } } ' | sed '/0\/1/d' | sed '/0\/2/d' | sed '/0\/3/d' | sed '/1\/2/d' | sed '/1\/3/d' | sed '/2\/3/d' > D4.snp.vcf



























# make the vcf into a table that includes only D species, not including D6, D10, D2-1, D2-2
# replace all unknown genotypes ("./.") with Q/Q so that we can parse in awk
# remove read counts and other metrics to leave only genotypes
# require that D4-Colima ($11) is different from D4 ($12) and that both have information
grep '^#CHROM\|^Chr' all.PASS_snp.vcf | awk ' BEGIN{OFS="\t";} { print $1,$2,$4,$5,$20,$21,$26,$30,$31,$32,$33,$34,$35,$36,$37,$38,$39,$40,$44,$45,$46,$47,$48  } ' | sed 's|[.]/[.]:|Q/Q:|g' |  awk -v OFS="\t"  ' { for ( i= 5; i<=23; i++) { gsub(":.*", "", $i) } print }  '  |  awk ' {  if ( ($11 != "Q/Q" && $12 != "Q/Q" && $11 != $12) || ($1 == "#CHROM")) { print } } ' > D4.snp.table

# if we want to test the D4s strictly against D3, the putative source of introgression
# remove all other D species
# require that all D3 species are uniform and homozygous at that SNP
awk ' { if ( (($8 == $9 && $8 == $10) && ($8 == "0/0" || $8 == "1/1" || $8 == "2/2" || $8 == "3/3")) || ($1 == "#CHROM" )) { print } } ' D4.snp.table | cut -f1,2,3,4,8,11,12 | sed 's/-D-27[.]Gdavi1//' | sed 's/[.]Garid1/C/' > D3D4.snp.table # 10,373,297


# comparing just D3, D4C, D4
# when D3=D4!=D4C -- this is the rate of autapomophy in D4C
# when D3=D4C!=D4 -- rate of autapomorphy in D4 + introgression
# the ratio of these should be 1:1 if there is no introgression

# first remove positions were D3 != D4C || D3 != D4
awk ' { if(( $5 == $6 ) || ( $5 == $7 ) || ( $1 == "#CHROM" )) { print } } ' D3D4.snp.table > D3vD4.snp.table # 9,579,600




# for positions homozygous in both D4 and D4C:
# if D3 == D4 != D4C; then rate of autapomorphy in D4C
# if D3 == D4C != D4; then homoplasious autapomorphy or autapomorphy in D4 or introgression
# all things considered equal, we expect the amount of D3 == D4 ~ D3 == D4C; divergence will suggest either increased autapomorphy in non-Colima D4 (unlikely to be substantial increase) or introgression
# have to consider heterozygosity in the != accession....



# for heterozygous positions in D4 or D4C
# cannot be heterozygous in both? but what is D3 == 0/0, D4C == 0/2, D4 == 1/2 -- what is our inference? shared autapomorphy in D4s plus introgression from D4C?
