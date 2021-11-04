plot.prop.mutation <- function(f,df){
	tmp <- as.data.frame.matrix(table(df[,f],df$Cluster))
	tmp$hLMS <- tmp$hLMS/28
	tmp$oLMS <- tmp$oLMS/11
	tmp$alteration <- row.names(tmp)
	tmp <- melt(tmp)
	tmp$value[tmp$value==0] <- NA
	tmp$alteration <- factor(as.vector(tmp$alteration),levels=c("2L","MS","MUT","MS/FS","MS/L","MUT/L","FS/L","L/WT","MUT/WT","FS/WT","WT","N","noyau","membrane","C","nothing","X"))
	colnames(tmp) <- c("alteration","Cluster","sample.freq")
	tmp$gene <- f
	tmp
}
