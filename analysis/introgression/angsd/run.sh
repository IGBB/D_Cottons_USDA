# curl -O https://raw.githubusercontent.com/ANGSD/angsd/master/R/jackKnife.R


angsd -doAbbababa 1 -doCounts 1 -anc D6.fa -uniqueOnly -minQ 20 -checkBamHeaders 0 -rf Chromosomes.txt -out aridum -bam bam.filelist
Rscript jackKnife.R file=aridum.abbababa indNames=bam.filelist outfile=aridum

parallel 'echo Chr{} > Chr{}.txt' ::: {01..13}
parallel 'angsd -doAbbababa 1 -doCounts 1 -anc D6.fa -uniqueOnly -minQ 20 -checkBamHeaders 0 -rf Chr{}.txt -out aridum{} -bam bam.filelist' ::: {01..13}
parallel 'Rscript jackKnife.R file=aridum{}.abbababa indNames=bam.filelist outfile=aridum{}' ::: {01..13}
