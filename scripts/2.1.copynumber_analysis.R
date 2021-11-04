message("Start Copy number analysis")
# Source necessary functions and libraries
source("2.functions.R")

chrom.size <- read.table(file.path(data.dir,"ext","chrom.sizes.hg19.tab"),sep="\t")
chrom.size$chrNum <- as.numeric(sub("chr","",chrom.size$V1))
chrom.size$chrNum[23:24] <- c(23,24)


load(file.path(data.dir,"RData","CN_merge_cohort.RData"))

load(file.path(data.dir,"RData","Matrice_Complexes_129_SNP6.Rdata"))
SNP6 <- dataG
load(file.path(data.dir,"RData","Matrice_Complexes_237_BAC.Rdata"))
BAC <- dataG

CGH <- merge(SNP6,BAC,by=0,all=T)
row.names(CGH) <- CGH$Row.names
CGH <- CGH[,-1]
CGH <- CGH[,colnames(CGH)%in%sub(".affy.sarcComp","",colnames(puces.expr))]

gene.symbols <- row.names(CGH)
CGH <- sapply(CGH,copy.test)
row.names(CGH) <- gene.symbols
na <- apply(CGH,1,function(x) length(which(is.na(x))))
CN.puces <- CGH[na<(ncol(CGH)/4),]

cytoband.hg19 <- import.cyto(file.path(data.dir,"ext"),Txdb,org.Hs.eg.db)

################################################################################@
# compare CN
message("Computing percentage events per cohort")
message("Affymetrix")
df.CN <- df.LMS[colnames(CN.puces),]
CN.genes.puces.noFDR <- import.copyNumber.by.genes(CN.puces,df.CN,cytoband.hg19,logRatio,0.01,fisher=FALSE)

message("TCGA")
df.CN.TCGA <- merge.clinical[row.names(merge.clinical)%in%colnames(CN_cohort) & merge.clinical$exp=="TCGA",]
CN_TCGA <- CN_cohort[,colnames(CN_cohort)%in%row.names(df.CN.TCGA)]
CN.genes.TCGA.noFDR <- import.copyNumber.by.genes(CN_TCGA,df.CN.TCGA,cytoband.hg19,logRatio.TCGA,0.01,fisher=FALSE)

message("ICGC")
df.CN.ICGC <- merge.clinical[row.names(merge.clinical)%in%colnames(CN_cohort) & merge.clinical$exp=="ICGC",]
CN_ICGC <- CN_cohort[,colnames(CN_cohort)%in%row.names(df.CN.ICGC)]
CN.genes.ICGC.noFDR <- import.copyNumber.by.genes(CN_ICGC,df.CN.ICGC,cytoband.hg19,logRatio.ICGC,0.01,fisher=FALSE)


message("Compare penetrance profiles between 3 cohorts per event types")
genes.dt.puces <- CN.genes.puces.noFDR$genes.dt
genes.dt.puces <- percentage.events(genes.dt.puces,cytoband.hg19)
m <- match(genes.dt.puces$gene,cytoband.hg19$SYMBOL)
genes.dt.puces$start.gene <- cytoband.hg19$gene.start[m]
genes.dt.puces <- genes.dt.puces[order(start.gene,decreasing=F)]
chrNum <- as.numeric(sub("chr","",as.vector(genes.dt.puces$chr)))
genes.dt.puces <- genes.dt.puces[order(chrNum,decreasing=F)]
genes.dt.puces$position <- 1:nrow(genes.dt.puces)

genes.dt.ICGC <- CN.genes.ICGC.noFDR$genes.dt
genes.dt.ICGC <- percentage.events(genes.dt.ICGC,cytoband.hg19)
chrNum <- as.numeric(sub("chr","",as.vector(genes.dt.ICGC$chr)))
m <- match(genes.dt.ICGC$gene,cytoband.hg19$SYMBOL)
genes.dt.ICGC$start.gene <- cytoband.hg19$gene.start[m]
genes.dt.ICGC <- genes.dt.ICGC[order(start.gene,decreasing=F)]
genes.dt.ICGC <- genes.dt.ICGC[order(chrNum,decreasing=F)]
genes.dt.ICGC$position <- 1:nrow(genes.dt.ICGC)

genes.dt.TCGA <- CN.genes.TCGA.noFDR$genes.dt
genes.dt.TCGA <- percentage.events(genes.dt.TCGA,cytoband.hg19)
m <- match(genes.dt.TCGA$gene,cytoband.hg19$SYMBOL)
genes.dt.TCGA$start.gene <- cytoband.hg19$gene.start[m]
genes.dt.TCGA <- genes.dt.TCGA[order(start.gene,decreasing=F)]
chrNum <- as.numeric(sub("chr","",as.vector(genes.dt.TCGA$chr)))
genes.dt.TCGA <- genes.dt.TCGA[order(chrNum,decreasing=F)]
genes.dt.TCGA$position <- 1:nrow(genes.dt.TCGA)

compare_3studies <- merge(genes.dt.TCGA[,list(gene,TCGA.gain.h=percent.gain.h,TCGA.lost.h=percent.lost.h,TCGA.amp.h=percent.amp.h,TCGA.del.h=percent.del.h,TCGA.gain.o=percent.gain.o,TCGA.lost.o=percent.lost.o,TCGA.amp.o=percent.amp.o,TCGA.del.o=percent.del.o)],
							genes.dt.ICGC[,list(gene,ICGC.gain.h=percent.gain.h,ICGC.lost.h=percent.lost.h,ICGC.amp.h=percent.amp.h,ICGC.del.h=percent.del.h,ICGC.gain.o=percent.gain.o,ICGC.lost.o=percent.lost.o,ICGC.amp.o=percent.amp.o,ICGC.del.o=percent.del.o)],by="gene",all=T)

compare_3studies <- merge(compare_3studies,
							genes.dt.puces[,list(gene,puces.gain.h=percent.gain.h,puces.lost.h=percent.lost.h,puces.amp.h=percent.amp.h,puces.del.h=percent.del.h,puces.gain.o=percent.gain.o,puces.lost.o=percent.lost.o,puces.amp.o=percent.amp.o,puces.del.o=percent.del.o)],by="gene",all=T)

Rsquared.all<- comp.rsquared(compare_3studies,"gain|lost")

Rsquared.all <- Rsquared.all[c(grep("gain.h",colnames(Rsquared.all)),grep("lost.h",colnames(Rsquared.all)),grep("gain.o",colnames(Rsquared.all)),grep("lost.o",colnames(Rsquared.all))),c(grep("gain.h",colnames(Rsquared.all)),grep("lost.h",colnames(Rsquared.all)),grep("gain.o",colnames(Rsquared.all)),grep("lost.o",colnames(Rsquared.all)))]
            
pdf(file.path(res.dir,"FS4A_heatmap_rsquared_CN.pdf"),width=5,heigh=4.7)
pheatmap(Rsquared.all,cluster_rows=F,cluster_cols=F,display_numbers=T,breaks=seq(-1,1,0.02),color=colorRampPalette(c(colours()[130],colours()[121],"white","white",colours()[404],"red"))(100))
dev.off()
#export.plot(file.path(res.dir,"FS4A_heatmap_rsquared_CN"),width=5,heigh=4.7)
message("Figure S4A OK")


################################################################################@
message("Computing percentage events and Fisher's exact test merged cohort")
df.CN_13 <- merge.clinical[row.names(merge.clinical)%in%colnames(CN_cohort),]
CN_13 <- CN_cohort[,na.omit(match(row.names(df.CN_13),colnames(CN_cohort)))]
df.CN_13$Cluster <- as.vector(df.CN_13$Cluster)

CN.all <- merge(CN_13,CN.puces,by=0,all=T)
row.names(CN.all) <- CN.all$Row.names
CN.all <- CN.all[,-1]
annots.all <- rbind(df.LMS[colnames(CN.puces),c("Location","Cluster")],df.CN_13[,c("Location","Cluster")])
logRatio.all <- merge(data.frame(TCGA=logRatio.TCGA,row.names=names(logRatio.TCGA)),data.frame(ICGC=logRatio.ICGC,row.names=names(logRatio.ICGC)),by=0,all=T)
m <- match(logRatio.all$Row.names,names(logRatio))
logRatio.all$puces <- logRatio[m]
logRatio.all <- data.table(logRatio.all)
logRatio.all$all <- apply(logRatio.all[,list(TCGA,ICGC,puces)],1,median,na.rm=T)
logRatio.mean <- logRatio.all$all
names(logRatio.mean) <- logRatio.all$Row.names

message("Computing percentage events and Fisher's exact test enrichment hLMS")
CN.genes.all.hLMS <- compute.CN(CN.all,annots.all,cytoband.hg19,logRatio.mean,logRatio,logRatio.TCGA,logRatio.ICGC,ref="hLMS")
genes.dt.all <- CN.genes.all.hLMS$genes.dt.all
fishers.band <- CN.genes.all.hLMS$fishers.band
annots.genes.all <- CN.genes.all.hLMS$annots.genes.all

message("Computing percentage events and Fisher's exact test enrichment oLMS")
CN.genes.all.oLMS <- compute.CN(CN.all,annots.all,cytoband.hg19,logRatio.mean,logRatio,logRatio.TCGA,logRatio.ICGC,ref="oLMS")
genes.dt.oLMS <- CN.genes.all.oLMS$genes.dt.all
fishers.band.oLMS <- CN.genes.all.oLMS$fishers.band

genes.dt.all$enr.hLMS <- fishers.band$enr[match(genes.dt.all$ont,fishers.band$ont)]
genes.dt.oLMS$enr.oLMS <- fishers.band.oLMS$enr[match(genes.dt.oLMS$ont,fishers.band.oLMS$ont)]

levels.cb <- as.vector(unique(genes.dt.all$ont))

fisher.events <- data.table(fisher.hLMS=genes.dt.all$fisher_genes,enr.hLMS=genes.dt.all$enr.hLMS,fisher.oLMS=genes.dt.oLMS$fisher_genes,enr.oLMS=genes.dt.oLMS$enr.oLMS,z.expr=genes.dt.all$z.expr,ont=genes.dt.all$ont)
fisher.events <- fisher.events[,list(nb.genes=length(fisher.hLMS),
									fisher.hLMS.pc=ifelse(unique(enr.hLMS)==0,0,sum(fisher.hLMS==unique(enr.hLMS))/length(fisher.hLMS)),
									fisher.oLMS.pc=ifelse(unique(enr.oLMS)==0,0,sum(fisher.oLMS==unique(enr.oLMS))/length(fisher.hLMS)),
									median.expr.diff= round(median(z.expr,na.rm=T),2),
									IQR=paste(round(quantile(z.expr,na.rm=T,probs=0.25),2),":",round(quantile(z.expr,na.rm=T,probs=0.75),2),sep=" ")),by=c("ont","enr.hLMS","enr.oLMS")]

pos.neg <- tolower(data.table(read.table(file.path(res.dir,"GSEA_res_positions_neg.tab"),sep="\t",head=T))[FDR<0.01]$Description)
pos.pos <- tolower(data.table(read.table(file.path(res.dir,"GSEA_res_positions_pos.tab"),sep="\t",head=T))[FDR<0.01]$Description)
fisher.events$pos.enr <- 0
fisher.events[ont%in%pos.neg]$pos.enr <- -1
fisher.events[ont%in%pos.pos]$pos.enr <- 1

percent.all <- rbindlist(apply(genes.dt.all,1,percentage.all))
genes.dt.all <- cbind(genes.dt.all,percent.all)

median.events <- genes.dt.all[,grep("percent|ont",colnames(genes.dt.all)),with=F][,lapply(.SD,median),by=ont]

band.events <- merge(fisher.events,median.events,by="ont")

tab.per.genes <- genes.dt.all[,list(cytogenic.band=ont,SYMBOL=gene,percent.del.h,percent.del.o,percent.del,percent.lost.h,percent.lost.o,percent.loss,percent.gain.h,percent.gain.o,percent.gain,percent.amp.h,percent.amp.o,percent.amp,fisher.h=fisher_genes,fisher.o=genes.dt.oLMS$fisher_genes,z.puces,z.TCGA,z.ICGC,z.median=z.expr)]
tab.band <- band.events[match(levels.cb,ont),list(cytogenic.band=ont,nb.genes,percent.del.h,percent.del.o,percent.del,percent.lost.h,percent.lost.o,percent.loss,percent.gain.h,percent.gain.o,percent.gain,percent.amp.h,percent.amp.o,percent.amp,enr.hLMS,fisher.hLMS.pc,enr.oLMS,fisher.oLMS.pc,GSEA=pos.enr,median.expr.diff,IQR)]

write.table(tab.band,file.path(res.dir,"TS3_CN_per_cytoband.tab"),sep="\t",quote=F,row.names=F)
write.table(tab.per.genes,file.path(res.dir,"TS3_CN_per_genes.tab"),sep="\t",quote=F,row.names=F)
message("TS3 OK")



################################################################################@
message("Visualize heatmap")
annoRow <- as.data.frame(annots.genes.all)
annoRow <- annoRow[annoRow$SYMBOL%in%genes.dt.all[CGH==TRUE]$gene,]
row.names(annoRow) <- annoRow$SYMBOL
annoRow <- annoRow[,c("fisher_genes","arm","chrNum")]
annoRow$chrNum <- factor(as.vector(annoRow$chrNum),levels=as.character(1:23))
colnames(annoRow)[1] <- "fisher.event"
CN <- CN.all[match(row.names(annoRow),row.names(CN.all)),]

gene.order <- order(as.numeric(as.vector(annoRow$chrNum)))

annoRow <- annoRow[,-3]


mat <- CN
chrNum.col <- rep(c("black","grey"),13)
names(chrNum.col) <- as.character(1:26)

            pdf(file.path(res.dir,"Figure3A_heatmap_CN_puce_available_genes.pdf"),width=6,heigh=6)
pheatmap(mat[gene.order,row.names(annots.all)[order(annots.all$Cluster)]],
	cluster_rows=F,cluster_cols=F,#res.hcpc$call$t$tree, 
	border_color=NA,
	color=colorRampPalette(c(colours()[131],"white","red"))(100), 
	show_rownames=F,show_colnames=F,breaks=seq(0,4,0.04),
	#gaps_row=chroms, 
	#gaps_col=cumsum(table(df.CN$Cluster)),
	annotation_row=annoRow,annotation_col=data.frame(Cluster=annots.all$Cluster,row.names=row.names(annots.all)),
	annotation_colors=c(generate.annot.col(annoRow[,c("fisher.event","arm")]),list(chrNum=chrNum.col,Cluster=c("hLMS"=colours()[613],"oLMS"=colours()[132]))),
	legend = FALSE)
dev.off()
#export.plot(file.path(res.dir,"Figure3A_heatmap_CN_puce_available_genes"),width=6,heigh=6)
message("Figure3A OK")

genes.dt.plot <- genes.dt.all[CGH==TRUE]
genes.dt.plot$position <- 1:nrow(genes.dt.plot)
chrNum.lines <- c(1,cumsum(Rle(genes.dt.plot$chrNum)@lengths))
chrNum.lines <- data.frame(x=chrNum.lines,xend=chrNum.lines,y=1,yend=-1)

message("Visualize CN hLMS penetrance")
g <- ggplot(genes.dt.plot,aes(x=position)) + 
	geom_hline(yintercept=0.5,linetype="dotted",size=0.1) +
	geom_hline(yintercept=-0.5,linetype="dotted",size=0.1) +
	geom_hline(yintercept=0,linetype="dashed",size=0.1) +
	geom_segment(aes(xend=position,y=0,yend=percent.gain.h),color=colours()[34],size=0.1) +
	geom_segment(aes(xend=position,y=percent.gain.h,yend=percent.gain.h+percent.amp.h),color=colours()[35],size=0.1) +
	geom_segment(aes(xend=position,y=0,yend=-percent.lost.h),color=colours()[121],size=0.1) +
	geom_segment(aes(xend=position,y=-percent.lost.h,yend=-percent.lost.h-percent.del.h),color=colours()[131],size=0.1) +
	#geom_segment(aes(x=x,xend=xend,y=y,yend=yend),data=chrNum.lines,size=0.1) +
	ylim(-1,1) + theme(panel.background = element_rect(fill = "white", colour = "grey")) 

            pdf(file.path(res.dir,"Figure3B_puces_avail_all_hLMS.pdf"),width=15,heigh=4)
print(g)
#export.plot(file.path(res.dir,"Figure3B_puces_avail_all_hLMS"),width=15,heigh=4)
dev.off()
message("Figure3B hLMS OK")

message("Visualize CN oLMS penetrance")
g <- ggplot(genes.dt.plot,aes(x=position)) + 
	geom_hline(yintercept=0.5,linetype="dotted",size=0.1) +
	geom_hline(yintercept=-0.5,linetype="dotted",size=0.1) +
	geom_hline(yintercept=0,linetype="dashed",size=0.1) +
	geom_segment(aes(xend=position,y=0,yend=percent.gain.o),color=colours()[34],size=0.1) +
	geom_segment(aes(xend=position,y=percent.gain.o,yend=percent.gain.o+percent.amp.o),color=colours()[35],size=0.1) +
	geom_segment(aes(xend=position,y=0,yend=-percent.lost.o),color=colours()[121],size=0.1) +
	geom_segment(aes(xend=position,y=-percent.lost.o,yend=-percent.lost.o-percent.del.o),color=colours()[131],size=0.1) +
	#geom_segment(aes(x=x,xend=xend,y=y,yend=yend),data=chrNum.lines,size=0.1) +
	ylim(-1,1) + theme(panel.background = element_rect(fill = "white", colour = "grey"))

pdf(file.path(res.dir,"Figure3B_puces_avail_all_oLMS.pdf"),width=15,heigh=4)
print(g)
#export.plot(file.path(res.dir,"Figure3B_puces_avail_all_oLMS"),width=15,heigh=4)
dev.off()
message("Figure3B oLMS OK")

####################################################################
# MYOCD 
message("MYOCD copy number and expression")
CN.MYOCD <- CN.all["MYOCD",row.names(annots.all)]
MYOCD.invest <- as.data.frame(t(CN.MYOCD))
MYOCD.invest$Cluster <- annots.all$Cluster
MYOCD.invest$exp <- c(rep("Affy",84),rep("ICGC",39),rep("TCGA",62))
MYOCD.expr <- cbind(puces.expr["MYOCD",colnames(CN.puces)],ICGC.expr["MYOCD",row.names(annots.all)[85:(85+38)]],TCGA.expr["MYOCD",colnames(CN_TCGA)])
MYOCD.invest$expression <- as.numeric(MYOCD.expr)
MYOCD.invest <- data.table(MYOCD.invest)

MYOCD.stats <- MYOCD.invest[,list(median.gain=median(expression[MYOCD>2],na.rm=T),median=median(expression,na.rm=T),Q1=quantile(expression,0.01,na.rm=T)),by=c("Cluster","exp")]

g <- ggplot(MYOCD.invest[MYOCD>2],aes(x=Cluster,y=expression)) + 
	facet_wrap(~exp,scales="free") + geom_boxplot(aes(color=Cluster)) + 
	theme_bw() + scale_color_manual(values=c(hLMS=colours()[613],oLMS=colours()[131],LMS="black")) +
	stat_compare_means(comparisons=list(c("hLMS","oLMS")))

pdf(file.path(res.dir,"FS4B_MYOCD_CN_amplified_expression.pdf"),width=6,heigh=6)
print(g)
#export.plot(file.path(res.dir,"FS4B_MYOCD_CN_amplified_expression"),width=6,heigh=6)
dev.off()
message("FS4B OK OK")



message("MYOCD expression ranking with inverse cumulative function")
MYOCD.invest <- as.data.frame(t(cbind(ICGC.expr["MYOCD",],TCGA.expr["MYOCD",],puces.expr["MYOCD",])))
cluster <- c(merge.clinical$Cluster,df.LMS$Cluster)
names(cluster) <- c(row.names(merge.clinical),sub("-",".",row.names(df.LMS)))
MYOCD.invest$cluster <- cluster[sub("[.]1$","",row.names(MYOCD.invest))]
MYOCD.invest$exp <- c(rep("ICGC",39),rep("TCGA",63),rep("Affy",98))

MYOCD.invest <- data.table(MYOCD.invest)

thr.affy <- quantile(as.numeric(as.matrix(puces.expr)),0.75)
thr.ICGC <- quantile(as.numeric(as.matrix(ICGC.expr)),0.75,na.rm=T)
thr.TCGA <- quantile(as.numeric(as.matrix(TCGA.expr)),0.75)

MYOCD.invest$MYOCD.quart <- "no"
MYOCD.invest[exp=="Affy" & MYOCD>=thr.affy]$MYOCD.quart <- "yes"
MYOCD.invest[exp=="TCGA" & MYOCD>=thr.TCGA]$MYOCD.quart <- "yes"
MYOCD.invest[exp=="ICGC" & MYOCD>=thr.ICGC]$MYOCD.quart <- "yes"

cumul.TCGA <- MYOCD.invest[exp=="TCGA"]
cumul.TCGA <- cumul.TCGA[order(MYOCD)]
cumul.TCGA$hLMS <- 1 - (cumsum(cumul.TCGA$cluster=="hLMS")/sum(cumul.TCGA$cluster=="hLMS"))
cumul.TCGA$oLMS <- 1 - (cumsum(cumul.TCGA$cluster=="oLMS")/sum(cumul.TCGA$cluster=="oLMS"))

pc.high.TCGA <- sum(cumul.TCGA[cluster=="hLMS"]$MYOCD.quart=="yes")/sum(cumul.TCGA$cluster=="hLMS")

g.TCGA <- ggplot(cumul.TCGA,aes(x=MYOCD,y=hLMS)) + geom_line(colour=colours()[613]) + 
	geom_line(aes(y=oLMS),colour=colours()[131]) + ylab("") + xlab("") +
	geom_hline(yintercept=pc.high.TCGA,size=0.3,linetype="dashed") +
	geom_vline(xintercept=thr.TCGA,size=0.3,linetype="dashed") + theme_bw() 


cumul.ICGC <- MYOCD.invest[exp=="ICGC"]
cumul.ICGC <- cumul.ICGC[order(MYOCD)]
cumul.ICGC$hLMS <- 1 - (cumsum(cumul.ICGC$cluster=="hLMS")/sum(cumul.ICGC$cluster=="hLMS"))
cumul.ICGC$oLMS <- 1 - (cumsum(cumul.ICGC$cluster=="oLMS")/sum(cumul.ICGC$cluster=="oLMS"))

pc.high.ICGC <- sum(cumul.ICGC[cluster=="hLMS"]$MYOCD.quart=="yes")/sum(cumul.ICGC$cluster=="hLMS")

g.ICGC <- ggplot(cumul.ICGC,aes(x=MYOCD,y=hLMS)) + geom_line(colour=colours()[613]) + 
	geom_line(aes(y=oLMS),colour=colours()[131]) + ylab("") + xlab("") +
	geom_hline(yintercept=pc.high.ICGC,size=0.3,linetype="dashed") +
	geom_vline(xintercept=thr.ICGC,size=0.3,linetype="dashed") + theme_bw()

cumul.Affy <- MYOCD.invest[exp=="Affy"]
cumul.Affy <- cumul.Affy[order(MYOCD)]
cumul.Affy$hLMS <- 1 - (cumsum(cumul.Affy$cluster=="hLMS")/sum(cumul.Affy$cluster=="hLMS"))
cumul.Affy$oLMS <- 1 - (cumsum(cumul.Affy$cluster=="oLMS")/sum(cumul.Affy$cluster=="oLMS"))

pc.high.Affy <- sum(cumul.Affy[cluster=="hLMS"]$MYOCD.quart=="yes")/sum(cumul.Affy$cluster=="hLMS")

g.Affy <- ggplot(cumul.Affy,aes(x=MYOCD,y=hLMS)) + geom_line(colour=colours()[613]) + 
	geom_line(aes(y=oLMS),colour=colours()[131]) + ylab("") + xlab("") +
	geom_hline(yintercept=pc.high.Affy,size=0.3,linetype="dashed") +
	geom_vline(xintercept=thr.affy,size=0.3,linetype="dashed") + theme_bw()

pdf(file.path(res.dir,"FS4C_cumul_MYOCD_expr.pdf"),width=6.5,heigh=2)
multiplot(g.Affy,g.ICGC,g.TCGA,cols=3)
#export.plot(file.path(res.dir,"FS4C_cumul_MYOCD_expr"),width=6.5,heigh=2)
dev.off()
message("FS4C OK")


message("Visualize zoom MYOCD hLMS penetrance")

m <- match(genes.dt.all$gene,cytoband.hg19$SYMBOL)
genes.dt.all$start.gene <- cytoband.hg19$gene.start[m]
#diff.genes <- row.names(df_tests_diff)[df_tests_diff$ttest.adj>2]
# 7452375 - 26658152
zoom.MYOCD <- genes.dt.all[13403:13561,list(SYMBOL=gene,ont,position=position,t.score.median=z.expr/max(abs(z.expr),na.rm=T),amp=percent.amp.h,gain=percent.gain.h,loss=-percent.lost.h,del=-percent.del.h)]
bands <- c(1,cumsum(Rle(zoom.MYOCD$ont)@lengths)) + 13402
zoom.MYOCD <- melt(zoom.MYOCD,id.vars=c("SYMBOL","position","ont"))
zoom.MYOCD$group <- ifelse(grepl("t.score.median",zoom.MYOCD$variable),"expression","penetrance")
g <- ggplot(zoom.MYOCD[zoom.MYOCD$group!="expression",],aes(x=position,y=value)) + geom_hline(yintercept=0,linetype="dashed",size=0.3) + 
	geom_col(aes(fill=variable),color="black",alpha=0.8) + 
	theme(panel.background = element_rect(fill = "white", colour = "black")) +#,panel.grid.major.y=element_line(size=0.1,linetype="dashed",colour="black")) +
	scale_fill_manual(values=c(loss=colours()[121],gain=colours()[34],amp=colours()[36],del=colours()[131])) +
	geom_segment(aes(xend=position,yend=0),data=zoom.MYOCD[zoom.MYOCD$group=="expression",],size=0.1) +
	geom_point(data=zoom.MYOCD[zoom.MYOCD$group=="expression",],aes(color=value)) + 
	scale_colour_gradient2(high="red",low="blue",mid="white",midpoint=0) + #geom_point(data=zoom.MYOCD[zoom.MYOCD$group=="expression" & zoom.MYOCD$SYMBOL%in%diff.genes,],shape=21,fill=NA,colour="black") +
  scale_y_continuous("value", sec.axis = sec_axis(~ . * 12.47618, name = "median t-score")) +
  	geom_segment(aes(x=band,xend=band,y=-1,yend=1),linetype="dashed",size=0.1,data=data.frame(band=bands))
	#geom_text_repel(aes(label=SYMBOL),data=zoom.MYOCD[zoom.MYOCD$group=="expression" & zoom.MYOCD$SYMBOL%in%diff.genes,])

pdf(file.path(res.dir,"Figure4A_zoom_MYOCD.pdf"),width=17,heigh=7)
print(g)
#export.plot(file.path(res.dir,"Figure4A_zoom_MYOCD"),width=17,heigh=7)
dev.off()
message("Figure 4A OK")

message("Compute MYOCD expression distribution")
MYOCD.expression <- rbind(data.frame(expression=as.numeric(puces.expr["MYOCD",]),Cluster=df.LMS$Cluster[na.omit(match(colnames(puces.expr),row.names(df.LMS)))],cohort="Affy"),
							data.frame(expression=as.numeric(ICGC.expr["MYOCD",colnames(ICGC.expr)%in%row.names(merge.clinical)]),Cluster=merge.clinical$Cluster[na.omit(match(colnames(ICGC.expr),row.names(merge.clinical)))],cohort="ICGC"),
							data.frame(expression=as.numeric(TCGA.expr["MYOCD",colnames(TCGA.expr)%in%row.names(merge.clinical)]),Cluster=merge.clinical$Cluster[na.omit(match(colnames(TCGA.expr),row.names(merge.clinical)))],cohort="TCGA"))

g <- ggplot(MYOCD.expression,aes(x=Cluster,y=expression)) + geom_violin(aes(color=Cluster)) +
	theme_pubr() + scale_color_manual(values=c(hLMS=colours()[613],oLMS=colours()[131])) + geom_jitter(width=0.1,aes(color=Cluster)) +
	facet_wrap(~cohort,scales="free") + stat_compare_means(comparison=list(c("hLMS","oLMS")),method="t.test")
pdf(file.path(res.dir,"Figure4B_MYOCD_expression_per_cohort.pdf"),width=5,heigh=3)
print(g)
#export.plot(file.path(res.dir,"Figure4B_MYOCD_expression_per_cohort"),width=5,heigh=3)
dev.off()
message("Figure 4B OK")



