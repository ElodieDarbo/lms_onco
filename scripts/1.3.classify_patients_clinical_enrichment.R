message("Start patient classification and clinical enrichment")
load(file.path(data.dir,"RData","merged_TCGA_ICGC_clinical.RData"))

ICGC.expr <- read.table(file.path(data.dir,"ext","ICGC_Matrice_Expression_FPKM_log2_2017_01_30.txt"),sep="\t",head=T,row=1)

TCGA.expr <- read.table(file.path(data.dir,"ext","TCGA_SARC_exp_HiSeqV2-2015-02-24.tab"),sep="\t",head=T,row=1)

message("Computing centroid from Affymetrix cohort")
common.genes <- row.names(puces.expr)[row.names(puces.expr)%in%row.names(na.omit(ICGC.expr))]
common.genes <- common.genes[common.genes%in%row.names(na.omit(TCGA.expr))]
puces.expr <- puces.expr[row.names(puces.expr)%in%common.genes,]
puces.diff <- puces.expr[row.names(puces.expr)%in%row.names(df_tests_diff)[df_tests_diff$ttest.adj>=2],]

centroid.LMS <- computeCentroids(df.LMS,puces.diff,"Cluster") 

puces.classif <- classify_LMS(puces.diff,centroid.LMS)
puces.clin <- df.LMS
puces.clin$centroid <- puces.classif$pred[match(row.names(puces.clin),row.names(puces.classif))]
puces.clin$hLMS.centroid <- puces.classif$hLMS[match(row.names(puces.clin),row.names(puces.classif))]
puces.clin$oLMS.centroid <- puces.classif$oLMS[match(row.names(puces.clin),row.names(puces.classif))]


message("Classify ICGC")
ICGC.classif <- classify_LMS(ICGC.expr,centroid.LMS)
ICGC.clin <- merge.clinical[merge.clinical$exp=="ICGC",]
ICGC.clin$centroid <- ICGC.classif$pred[match(row.names(ICGC.clin),row.names(ICGC.classif))]
ICGC.clin$hLMS.centroid <- ICGC.classif$hLMS[match(row.names(ICGC.clin),row.names(ICGC.classif))]
ICGC.clin$oLMS.centroid <- ICGC.classif$oLMS[match(row.names(ICGC.clin),row.names(ICGC.classif))]
ICGC.clin <- ICGC.clin[!is.na(ICGC.clin$centroid),]

message("Classify TCGA")
TCGA.classif <- classify_LMS(TCGA.expr,centroid.LMS)
TCGA.clin <- merge.clinical[merge.clinical$exp=="TCGA",]
TCGA.clin$centroid <- TCGA.classif$pred[match(row.names(TCGA.clin),row.names(TCGA.classif))]
TCGA.clin$hLMS.centroid <- TCGA.classif$hLMS[match(row.names(TCGA.clin),row.names(TCGA.classif))]
TCGA.clin$oLMS.centroid <- TCGA.classif$oLMS[match(row.names(TCGA.clin),row.names(TCGA.classif))]
TCGA.clin <- TCGA.clin[!is.na(TCGA.clin$centroid),]

merge.all.clin <- rbind(ICGC.clin[,c("hLMS.centroid","oLMS.centroid","centroid")],TCGA.clin[,c("hLMS.centroid","oLMS.centroid","centroid")],puces.clin[,c("hLMS.centroid","oLMS.centroid","centroid")])

message("Compute mclust gaussian mixture")
bimod.all <- merge.all.clin$hLMS.centroid
bimod.all.gmm = Mclust(bimod.all)
summary(bimod.all.gmm)

pdf(file.path(res.dir,"Figure1F_centroid_classification_mclust.pdf"),width=8,heigh=8)
par(mfrow=c(2,2))
plot(bimod.all.gmm,what="BIC")
plot(bimod.all.gmm,what="classification")
plot(bimod.all.gmm,what="uncertainty")
plot(bimod.all.gmm,what="density")
par(mfrow=c(1,1))
dev.off()
#export.plot(file.path(res.dir,"Figure1F_centroid_classification_mclust"),width=8,heigh=8)
message("Figure1F OK")

merge.all.clin$bimod.classif <- bimod.all.gmm$classification

merge.all.clin$classif.bimod <- "unclass"
merge.all.clin[merge.all.clin$bimod.classif==1,]$classif.bimod <- "hLMS"
merge.all.clin[merge.all.clin$bimod.classif==3,]$classif.bimod <- "oLMS"

range(merge.all.clin[merge.all.clin$classif.bimod=="hLMS",]$hLMS.centroid)
#0.1190197 0.6261260
range(merge.all.clin[merge.all.clin$classif.bimod=="unclass",]$hLMS.centroid)
#0.6390073 1.3872967
range(merge.all.clin[merge.all.clin$classif.bimod=="oLMS",]$hLMS.centroid)
#1.438048 1.870701

m <- match(row.names(merge.clinical),row.names(merge.all.clin))
merge.clinical$Cluster <- merge.all.clin$classif.bimod[m]
merge.clinical$hLMS.centroid <- merge.all.clin$hLMS.centroid[m]
merge.clinical$oLMS.centroid <- merge.all.clin$oLMS.centroid[m]
merge.clinical$centroid <- merge.all.clin$centroid[m]
merge.clinical$classif.bimod <- merge.all.clin$classif.bimod[m]
merge.clinical <- merge.clinical[merge.clinical$Cluster!="unclass",]

###############################################################################
message("Compute centroid for cell lines")

OC80 <- read.table(file.path(data.dir,"precomp","OC80_FPKM_log2.txt"),sep="\t",row=1)
OC48 <- read.table(file.path(data.dir,"precomp","OC48_FPKM_log2.txt"),sep="\t",row=1)
OC88 <- read.table(file.path(data.dir,"precomp","OC88_FPKM_log2.txt"),sep="\t",row=1)
m <- match(row.names(puces.diff),row.names(OC80))
OC80 <- data.table(SYMBOL=row.names(OC80)[na.omit(m)],OC80=OC80[na.omit(m),])
m <- match(row.names(puces.diff),row.names(OC48))
OC48 <- data.table(SYMBOL=row.names(OC48)[na.omit(m)],OC48=OC48[na.omit(m),])
m <- match(row.names(puces.diff),row.names(OC88))
OC88 <- data.table(SYMBOL=row.names(OC88)[na.omit(m)],OC88=OC88[na.omit(m),])
OCs <- data.frame(OC80=OC80$OC80,OC48=OC48$OC48,OC88=OC88$OC88,row.names=OC80$SYMBOL)
OCs.classif <- classify_LMS(OCs,centroid.LMS)
OCs.classif$OC <- row.names(OCs.classif)
OCs.classif$distance <- ifelse(OCs.classif$pred=="hLMS",OCs.classif$hLMS,-(1-OCs.classif$oLMS))
OCs.classif$OC <- factor(as.vector(OCs.classif$OC),levels=c("OC80","OC88","OC48"))
ggplot(OCs.classif,aes(x=OC,y=distance)) + geom_col(aes(fill=pred)) + theme_light() + ylim(-1,1) +
	scale_fill_manual(values=c(oLMS=colours()[121],hLMS=colours()[613],unclass="grey"))

pdf(file.path(res.dir,"Figure5A_OCs_centroids.pdf"),width=4,heigh=5)
print(ggplot(OCs.classif,aes(x=OC,y=distance)) + geom_col(aes(fill=pred)) + theme_light() + ylim(-1,1) +
	scale_fill_manual(values=c(oLMS=colours()[121],hLMS=colours()[613]))
)
dev.off()
###############################################################################


message("Compute Survival")
data.surv <- df.LMS[!is.na(df.LMS$Time_MFS),]
data.surv$Cluster <- factor(data.surv$Cluster)
data.surv$Meta <- as.numeric(as.character(data.surv$Meta))
fit <- survfit(Surv(Time_MFS,Meta) ~ Cluster,data=data.surv)
pdf(file.path(res.dir,"Figure1C_Affy_MFS_surv.pdf"),width=4,heigh=6)
print(ggsurvplot(fit, data = data.surv,  size = 1,  pval = TRUE,palette=c(colours()[613],colours()[131]),risk.table = TRUE,  risk.table.col = "strata", risk.table.height = 0.25))
#export.plot(file.path(res.dir,"Figure1C_Affy_MFS_surv"),width=4,heigh=6)
dev.off()
message("Figure1C OK")

data.surv <- merge.clinical[merge.clinical$exp=="ICGC",]
data.surv$Meta <- as.numeric(as.character(data.surv$Meta))
fit <- survfit(Surv(Time_MFS,Meta) ~ Cluster,data=data.surv)
print(ggsurvplot(fit, data = data.surv,  size = 1,  pval = TRUE,palette=c(colours()[613],colours()[131]),risk.table = TRUE,  risk.table.col = "strata", risk.table.height = 0.25))
#export.plot(file.path(res.dir,"ICGC_MSF_surv"))

data.surv <- merge.clinical[merge.clinical$exp=="TCGA",]
data.surv$Meta <- as.numeric(as.character(data.surv$Meta))
fit <- survfit(Surv(Time_MFS,Meta) ~ Cluster,data=data.surv)
print(ggsurvplot(fit, data = data.surv,  size = 1,  pval = TRUE,palette=c(colours()[613],colours()[131]),risk.table = TRUE,  risk.table.col = "strata", risk.table.height = 0.25))
#export.plot(file.path(res.dir,"TCGA_MFS_surv"))


#####################################################
message("Compute Clinical enrichment")

colnames(df.LMS)[colnames(df.LMS)=="Localisation"] <- "Location"
features <- c("differentiation","Location","grade_primary_tumour","sex")

fisher_clinical <- rbindlist(lapply(features,clinical.fisher,df.LMS=df.LMS))
cl <- df.LMS$Cluster[!is.na(df.LMS$mitotic_count_10_HPF)]
mitotic_count_10_HPF <- df.LMS$mitotic_count_10_HPF[!is.na(df.LMS$mitotic_count_10_HPF)]
wilcox_mitotic <- wilcox.test(mitotic_count_10_HPF[cl=="hLMS"],mitotic_count_10_HPF[cl=="oLMS"])

tmp <- data.frame(feat1="median mitotic counts cluster1",feat2="median mitotic counts cluster2",
	testfisher=wilcox_mitotic$p.value,
	nb.feat1.cl1=median(mitotic_count_10_HPF[cl=="hLMS"]),
	nb.feat1.cl2=median(mitotic_count_10_HPF[cl=="oLMS"]),
	nb.feat2.cl1=median(mitotic_count_10_HPF[cl=="oLMS"]),
	nb.feat2.cl2=median(mitotic_count_10_HPF[cl=="hLMS"]),
	feature="Mitotic counts")

fisher_clinical <- rbind(fisher_clinical,tmp)
write.table(fisher_clinical,file.path(res.dir,"Table1_Affymetrix_clinical.tab"),sep="\t",row.names=F,quote=F)

message("Table 1 OK")


#####################################################
features <- c("Location","Sex","Differentiation")
fisher_clinical <- rbindlist(lapply(features,clinical.fisher,merge.clinical))

cl <- merge.clinical$Cluster[!is.na(merge.clinical$Mitotic_count_10_HPF)]
Mitotic_count_10_HPF <- merge.clinical$Mitotic_count_10_HPF[!is.na(merge.clinical$Mitotic_count_10_HPF)]
wilcox_mitotic <- wilcox.test(Mitotic_count_10_HPF[cl=="hLMS"],Mitotic_count_10_HPF[cl=="oLMS"])


tmp <- data.frame(feat1="median mitotic counts cluster2",
	feat2="median mitotic counts cluster3",
	testfisher=wilcox_mitotic$p.value,
	nb.feat1.cl1=median(mitotic_count_10_HPF[cl=="hLMS"]),
	nb.feat1.cl2=median(mitotic_count_10_HPF[cl=="oLMS"]),
	nb.feat2.cl1=median(mitotic_count_10_HPF[cl=="oLMS"]),
	nb.feat2.cl2=median(mitotic_count_10_HPF[cl=="hLMS"]),
	feature="Mitotic counts")

fisher_clinical <- rbind(fisher_clinical,tmp)
write.table(fisher_clinical,file.path(res.dir,"Table2_ICGC_TCGA_clinical.tab"),sep="\t",row.names=F,quote=F)
message("Table 2 OK")

#####################################################
message("Compute gene expression variance")

sig.genes.affy <- row.names(df_tests_diff)[df_tests_diff$ttest.adj>=2]

var.hLMS.affy <- apply(puces.expr[row.names(puces.expr)%in%sig.genes.affy,row.names(df.LMS)[df.LMS$Cluster=="hLMS"]],1,var)
var.oLMS.affy <- apply(puces.expr[row.names(puces.expr)%in%sig.genes.affy,row.names(df.LMS)[df.LMS$Cluster=="oLMS"]],1,var)

var.hLMS.ICGC <- apply(ICGC.expr[row.names(ICGC.expr)%in%sig.genes.affy,colnames(ICGC.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="hLMS"]],1,var)
var.oLMS.ICGC <- apply(ICGC.expr[row.names(ICGC.expr)%in%sig.genes.affy,colnames(ICGC.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="oLMS"]],1,var)

var.hLMS.TCGA <- apply(TCGA.expr[row.names(TCGA.expr)%in%sig.genes.affy,colnames(TCGA.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="hLMS"]],1,var)
var.oLMS.TCGA <- apply(TCGA.expr[row.names(TCGA.expr)%in%sig.genes.affy,colnames(TCGA.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="oLMS"]],1,var)

var.tab <- rbind(data.table(variance=var.hLMS.affy,Cluster="hLMS Affy",group="hLMS",cohort="Affymetrix"),
				data.table(variance=var.oLMS.affy,Cluster="oLMS Affy",group="oLMS",cohort="Affymetrix"),
				data.table(variance=var.hLMS.ICGC,Cluster="hLMS ICGC",group="hLMS",cohort="ICGC"),
				data.table(variance=var.oLMS.ICGC,Cluster="oLMS ICGC",group="oLMS",cohort="ICGC"),
				data.table(variance=var.hLMS.TCGA,Cluster="hLMS TCGA",group="hLMS",cohort="TCGA"),
				data.table(variance=var.oLMS.TCGA,Cluster="oLMS TCGA",group="oLMS",cohort="TCGA"))

g <- ggplot(var.tab,aes(x=group,y=variance)) + geom_boxplot(aes(colour=group)) + scale_colour_manual(values=c("hLMS"=colours()[613],oLMS=colours()[132],gLMS=colours()[121])) +
	stat_compare_means(comparison=list(c("hLMS","oLMS")),method="wilcox") +
	facet_wrap(~cohort,scales="free") + theme(text=element_text(size=16)) + theme_bw()

pdf(file.path(res.dir,"boxplot_variances_all_cohorts_wilcox.pdf"),width=6,heigh=4)
print(g)
dev.off()
#export.plot(file.path(res.dir,"boxplot_variances_all_cohorts_wilcox"),width=6,height=3)
message("boxplot from variance per cohort")

#####################################################
message("Compute correlation between t-scores 3 cohorts")

logRatio <- df_tests_diff$logRatio
names(logRatio) <- row.names(df_tests_diff)

TCGA.expr <- TCGA.expr[,na.omit(match(row.names(merge.clinical[merge.clinical$exp=="TCGA",]),colnames(TCGA.expr)))]
ttest.TCGA <- rbindlist(apply(TCGA.expr,1,t.score,TCGA.annots=merge.clinical[merge.clinical$exp=="TCGA",]))
ttest.TCGA$SYMBOL <- row.names(TCGA.expr)
ttest.TCGA$ttest.adj.TCGA <- -log10(p.adjust(ttest.TCGA$ttest))
logRatio.TCGA <- ttest.TCGA$tscore
names(logRatio.TCGA) <- row.names(TCGA.expr)

ICGC.expr <- ICGC.expr[,na.omit(match(row.names(merge.clinical[merge.clinical$exp=="ICGC",]),colnames(ICGC.expr)))]
ttest.ICGC <- rbindlist(apply(ICGC.expr,1,t.score,TCGA.annots=merge.clinical[merge.clinical$exp=="ICGC",]))
ttest.ICGC$SYMBOL <- row.names(ICGC.expr)
ttest.ICGC$ttest.adj.ICGC <- -log10(p.adjust(ttest.ICGC$ttest))
logRatio.ICGC <- ttest.ICGC$tscore
names(logRatio.ICGC) <- row.names(ICGC.expr)

m1 <- match(names(logRatio),names(logRatio.TCGA))
m2 <- match(names(logRatio),names(logRatio.ICGC))
tscores <- data.frame(Affy=logRatio,ICGC=logRatio.ICGC[m2],TCGA=logRatio.TCGA[m1])
tscores$signature <- ifelse(row.names(tscores)%in%sig.genes.affy,"yes","no")

fit.AI <- round(summary(lm(Affy~ICGC,tscores))$r.squared,2)
fit.AT <- round(summary(lm(Affy~TCGA,tscores))$r.squared,2)
fit.IT <- round(summary(lm(ICGC~TCGA,tscores))$r.squared,2)

gAI <- ggplot(tscores[order(tscores$signature),],aes(x=Affy,y=ICGC)) + 
	geom_point(aes(color=signature),size=0.1) + 
	geom_smooth(method="lm",size=0.4) +
	scale_color_manual(values=c(yes="red",no="black")) +
	theme_bw() + ggtitle(paste("Affymetrix vs ICGC: R2=",fit.AI))
gAT <- ggplot(tscores[order(tscores$signature),],aes(x=Affy,y=TCGA)) + 
	geom_point(aes(color=signature),size=0.1) + 
	geom_smooth(method="lm",size=0.4) +
	scale_color_manual(values=c(yes="red",no="black"))+
	theme_bw() + ggtitle(paste("Affymetrix vs TCGA: R2=",fit.AT))
gIT <- ggplot(tscores[order(tscores$signature),],aes(x=ICGC,y=TCGA)) + 
	geom_point(aes(color=signature),size=0.1) + 
	geom_smooth(method="lm",size=0.4) +
	scale_color_manual(values=c(yes="red",no="black")) +
	theme_bw() + ggtitle(paste("ICGC vs TCGA: R2=",fit.IT))

pdf(file.path(res.dir,"tScores_3_cohorts_comparison.pdf"),width=12,heigh=4)
multiplot(gAI,gAT,gIT,cols=3)
#export.plot(file.path(res.dir,"tScores_3_cohorts_comparison"),width=12,heigh=4)
dev.off()
message("scatter plot from variance pairwise comparison")

var.hLMS.ICGC <- apply(ICGC.expr[,colnames(ICGC.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="hLMS"]],1,var)
var.oLMS.ICGC <- apply(ICGC.expr[,colnames(ICGC.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="oLMS"]],1,var)

var.hLMS.TCGA <- apply(TCGA.expr[,colnames(TCGA.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="hLMS"]],1,var)
var.oLMS.TCGA <- apply(TCGA.expr[,colnames(TCGA.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="oLMS"]],1,var)

median.hLMS.ICGC <- apply(ICGC.expr[,colnames(ICGC.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="hLMS"]],1,median)
median.oLMS.ICGC <- apply(ICGC.expr[,colnames(ICGC.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="oLMS"]],1,median)

median.hLMS.TCGA <- apply(TCGA.expr[,colnames(TCGA.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="hLMS"]],1,median)
median.oLMS.TCGA <- apply(TCGA.expr[,colnames(TCGA.expr)%in%row.names(merge.clinical)[merge.clinical$Cluster=="oLMS"]],1,median)

all.equal(colnames(puces.expr),row.names(df.LMS))
median.hLMS.puces <- apply(puces.expr[,colnames(puces.expr)%in%row.names(df.LMS)[df.LMS$Cluster=="hLMS"]],1,median)
median.oLMS.puces <- apply(puces.expr[,colnames(puces.expr)%in%row.names(df.LMS)[df.LMS$Cluster=="oLMS"]],1,median)


stats.ICGC <- data.frame(SYMBOL=names(var.oLMS.ICGC),var.hLMS.ICGC,var.oLMS.ICGC,median.hLMS.ICGC,median.oLMS.ICGC)
stats.ICGC <- merge(ttest.ICGC[,list(SYMBOL,ttest.adj.ICGC,tscore.ICGC=tscore)],stats.ICGC,by="SYMBOL")

stats.TCGA <- data.frame(SYMBOL=names(var.oLMS.TCGA),var.hLMS.TCGA,var.oLMS.TCGA,median.hLMS.TCGA,median.oLMS.TCGA)
stats.TCGA <- merge(ttest.TCGA[,list(SYMBOL,ttest.adj.TCGA,tscore.TCGA=tscore)],stats.TCGA,by="SYMBOL")

stats.puces <- data.frame(SYMBOL=names(median.hLMS.puces),median.hLMS.puces,median.oLMS.puces)
ttest.puces <- data.frame(SYMBOL=row.names(df_tests_diff),ttest.adj.puces=df_tests_diff$ttest.adj,tscore.puces=df_tests_diff$logRatio,var.hLMS.puces=df_tests_diff$var.hLMS,var.oLMS.puces=df_tests_diff$var.oLMS)
stats.puces <- merge(ttest.puces,stats.puces,by="SYMBOL")

stats.puces$signature <- ifelse(stats.puces$ttest.adj.puces>=2,1,0)

stats.all <- merge(stats.puces,stats.ICGC,by="SYMBOL",all=T)
stats.all <- merge(stats.all,stats.TCGA,by="SYMBOL",all=T)
write.table(stats.all,file.path(res.dir,"TS2A_stat_table_all_cohorts.tab"),sep="\t",row.names=F,quote=F)
message("TS2A OK")








