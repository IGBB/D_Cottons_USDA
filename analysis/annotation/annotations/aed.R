
args = commandArgs(trailingOnly=TRUE)
data <- read.table(args[1])
name <- args[2]
library(ggplot2)
ggplot(data, aes(V1))+stat_ecdf(geom="line")
ggsave(paste(name,".aed.png", sep=""))
