df <- read.table("rates.all", header=F)
names(df) <- c("node", "loss", "dup", "gain")
list <- unique(df$node)
list <- unique(grep("_|base|H", df$node, value=T, invert=T))
df <- df[df$node %in% list,]
df$node <- factor(df$node)

vals <- c("std_loss","std_gain")
df[,vals] <- 0

real <- read.table("calculated.rates", header=T, sep="\t")
row.names(real) <- real$node

for (row in c(1:nrow(df))) {
    for (names in row.names(real)) {
        ifelse(df$node[row]==names, df$std_loss[row] <- df$loss[row]/real[names,"mya"], df$std_loss[row] <- df$std_loss[row]) 
        ifelse(df$node[row]==names, df$std_gain[row] <- df$gain[row]/real[names,"mya"], df$std_gain[row] <- df$std_gain[row]) 
	}
}



###### losses ######
boxplot(std_loss~node, data=df, ylim=c(0,0.3))

points(factor(real$node), real$std_loss, col=2)


vals <- c("min", "max", "mean", "within")
real[,vals] <- "na"
row.names(real) <- real$node
real$node <- NULL


for (node in (row.names(real))) {
    real[node,"min"] <- round(min(df[df$node==node,]$loss), digits=3)
    real[node,"max"] <- round(max(df[df$node==node,]$loss), digits=3)  
	real[node,"mean"] <- round(mean(df[df$node==node,]$loss), digits=3)
	ifelse(real[node,"loss"] < real[node,"min"], real[node,"within"] <- "too low", 
	    ifelse(real[node,"loss"] > real[node,"max"], real[node,"within"] <- "too high", 
		real[node,"within"] <- "good"))
}

###### gains ######
boxplot(std_gain~node, data=df)

points(factor(real$node), real$std_gain, col=2)

###### duplications ######
boxplot(std_dup~node, data=df)

points(factor(real$node), real$std_dup, col=2)