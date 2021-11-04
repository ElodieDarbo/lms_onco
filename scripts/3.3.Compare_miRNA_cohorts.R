# This script aims at compare miRNA results from ICGC and TCGA

message("Start comparing ICGC/TCGA miRNA expression")

annots.TCGA <- merge.clinical[merge.clinical$exp=="TCGA",]
annots.ICGC <- merge.clinical[merge.clinical$exp=="ICGC",]

## Figure 2A: logFC comparison
# merge data from both cohorts
miRNA_ICGC_TCGA <- merge(miRNA.ICGC[,list(ID,hLMS.ICGC=median.hLMS,oLMS.ICGC=median.oLMS,logFC.ICGC=logFC,adj.pval.ICGC=adj.pval)],
                         miRNA.diffexp[,list(ID,hLMS.TCGA=median.hLMS,oLMS.TCGA=median.oLMS,logFC.TCGA=logratio,adj.pval.TCGA=adj.pval)],by="ID")

# Keep miRNA expressed (median norm counts > 1) in at least one group in each cohort
miRNA_ICGC_TCGA <- miRNA_ICGC_TCGA[((hLMS.ICGC>1 | oLMS.ICGC>1) ) & ((hLMS.TCGA>1 | oLMS.TCGA>1))]

# Add information about DE significance
miRNA_ICGC_TCGA$Significance <- "no"
miRNA_ICGC_TCGA[adj.pval.TCGA<0.01 & abs(logFC.TCGA) > 1  ]$Significance <- "TCGA"
miRNA_ICGC_TCGA[adj.pval.ICGC<0.01 & abs(logFC.ICGC) > 1]$Significance <- "ICGC"
miRNA_ICGC_TCGA[adj.pval.TCGA<0.01 & adj.pval.ICGC<0.01 & abs(logFC.TCGA) > 1 & abs(logFC.ICGC) > 1]$Significance <- "TCGA+ICGC"
#save(miRNA_ICGC_TCGA,file=file.path(res.dir,"merged_TCGA_ICGC_DE_miRNA.RData"))
#write.table(miRNA_ICGC_TCGA,file.path(res.dir,"ICGC_TCGA_miRNA_DE.tab"),sep="\t",row.names=F,quote=F)
message("Plot logFC scatterplot")
g <- ggplot(miRNA_ICGC_TCGA[order(Significance,decreasing=T)],aes(x=logFC.TCGA,y=logFC.ICGC)) + 
  geom_point(aes(color=Significance),alpha=0.7) + 
  geom_hline(yintercept = 0,linetype="dashed",size=0.3) + 
  geom_vline(xintercept = 0,linetype="dashed",size=0.3) +
  geom_smooth(method="lm",color="black",size=0.4) +
  theme_bw(base_rect_size = 1) + 
  scale_colour_manual(values=c(no="grey",ICGC="grey",TCGA="grey","TCGA+ICGC"="red"))
pdf(file.path(res.dir,"Fgure2A_XY_logFC_ICGC_TCGA_miRNA.pdf"),width=5,heigh=4)
print(g)
#export.plot(file.path(res.dir,"Fgure2A_XY_logFC_ICGC_TCGA_miRNA"),width=5,heigh=4)
dev.off()
message("Figure2A OK")

# compute R-squared value
message("Computing linear regression R2 score")
fit <- lm(logFC.TCGA~logFC.ICGC,data=miRNA_ICGC_TCGA)
print(summary(fit))


# Look at DIO3-DLK1 miRNA cluster
message("Analyse expression of not differentially expressed DIO3-DLK1 miRNAs")
load(file.path(data.dir,"RData","DIO3_DLK1_miRNA.RData")) # DIO3

miRNA_ICGC_TCGA$hairpin <- sub("-[35]p","",as.vector(miRNA_ICGC_TCGA$ID))

m <- match(tolower(miRNA_ICGC_TCGA$hairpin),DIO3)
miRNA_ICGC_TCGA$DIO3_DLK1 <- ifelse(is.na(m),"no","yes")

### Figure 2C: DIO3-DLK1 miRNA no significant DE 
miRNA_ICGC_TCGA_melt <- melt(miRNA_ICGC_TCGA[DIO3_DLK1=="yes" & Significance!="TCGA+ICGC",list(hLMS.ICGC,oLMS.ICGC,hLMS.TCGA,oLMS.TCGA,logFC.TCGA,logFC.ICGC)])
colnames(miRNA_ICGC_TCGA_melt) <- c("group","value")
miRNA_ICGC_TCGA_melt <- data.table(miRNA_ICGC_TCGA_melt)
miRNA_ICGC_TCGA_melt$LMS <- ifelse(grepl("hLMS",miRNA_ICGC_TCGA_melt$group),"hLMS","oLMS")
miRNA_ICGC_TCGA_melt$LMS[grepl("logFC",miRNA_ICGC_TCGA_melt$group)] <- "logFC"
miRNA_ICGC_TCGA_melt$feature <- ifelse(grepl("logFC",miRNA_ICGC_TCGA_melt$group),"logFC","median")
miRNA_ICGC_TCGA_melt$cohort <- ifelse(grepl("TCGA",miRNA_ICGC_TCGA_melt$group),"TCGA","ICGC")

write.table(miRNA_ICGC_TCGA,file.path(res.dir,"miRNA_ICGC_TCGA_DE.tab"),sep="\t",row.names=F,quote=F)
message("Table S2B OK")

g <- ggplot(miRNA_ICGC_TCGA_melt,aes(x=cohort,y=value)) + geom_boxplot(aes(colour=LMS)) + 
  scale_colour_manual(values=c("hLMS"=colours()[613],"oLMS"=colours()[132],logFC="black")) +
  facet_wrap(~feature,ncol=2,scales="free") + theme_bw()
pdf(file.path(res.dir,"FigureS3AB_boxplots_DIO3_DLK1_no_sig.pdf"),width=8,heigh=4)
print(g)
#export.plot(file.path(res.dir,"FigureS3AB_boxplots_DIO3_DLK1_no_sig"),width=8,heigh=4)
dev.off()
message("Figure S3A (logFC) et S3B (median) OK")


