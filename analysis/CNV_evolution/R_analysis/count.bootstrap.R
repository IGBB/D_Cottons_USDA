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


###### violin plot ###### 


#df$node<- as.factor(df$node)
#
#a <- ggplot(df, aes(x=node, y=std_gain)) + geom_violin()+ geom_boxplot(width=0.1, fill="blue")+ ylim(0,0.35)+coord_flip()+geom_point(data=real, aes(x=node, y=std_loss), size=2)
#b <- ggplot(df, aes(x=node, y=std_gain)) + geom_violin()+ geom_boxplot(width=0.1, fill="blue")+ ylim(2,4.25)+coord_flip()+geom_point(data=real, aes(x=node, y=std_loss), size=2)
#c <- ggplot(df, aes(x=node, y=std_loss)) + geom_violin()+ geom_boxplot(width=0.1, fill="blue")+ ylim(0,1)+ coord_flip()+geom_point(data=real, aes(x=node, y=std_loss), size=2)
#d <- ggplot(df, aes(x=node, y=std_loss)) + geom_violin()+ geom_boxplot(width=0.1, fill="blue")+ ylim(1,6.5)+ coord_flip()+geom_point(data=real, aes(x=node, y=std_loss), size=2)

#png("gain.loss.violin.png", 10000, 7500, pointsize=12, res=600)
#
#grid.arrange(a,b,c,d,ncol=2)

#dev.off()




###### losses ######

png("loss.boxplot.png", 10000, 7500, pointsize=12, res=600)
boxplot(std_loss~node, data=df, ylim=c(0,6.1), boxwex=0.2)
points(factor(real$node), real$std_loss, col=51, pch=24)
dev.off()


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
png("gain.boxplot.png", 10000, 7500, pointsize=12, res=600)
boxplot(std_gain~node, data=df, ylim=c(0,4.2), boxwex=0.2)
points(factor(real$node), real$std_gain, col=51, pch=24)
dev.off()

###### duplications ######
boxplot(std_dup~node, data=df)

points(factor(real$node), real$std_dup, col=2)