
chitable <- read.table("shared.derived.snps", sep="\t", header=TRUE, row.names=1)
chitable$p.value <- apply(chitable, 1, function(x) chisq.test(x, simulate.p.value=TRUE)$p.value)
