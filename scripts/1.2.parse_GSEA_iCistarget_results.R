
message("Start functional enrichment visualisation")

# Source necessary functions and libraries
source("1.functions.R")

##############################################################################################
# Parsing enrichment results
# GSEA
load(file.path(res.dir,"DE_LMS_affy_objects.RData"))
message("GSEA")

gsea.dir <- file.path(data.dir,"precomp","GSEA")
GSEA.files <- system(paste("ls",file.path(gsea.dir,"*","gsea_report_for_na_*html")),intern=T)
GSEA_res <- rbindlist(lapply(GSEA.files,readGSEA,gsea.dir=gsea.dir,res.dir=res.dir,prefix="GSEA_res"))

GSEA_res <- GSEA_res[FDR<0.01]
GSEA_res$signal <-  as.numeric(sub("%","",unlist(sapply(strsplit(as.vector(GSEA_res$signal),", signal="),"[[",2))))
GSEA_res$NES <- as.numeric(as.vector(GSEA_res$NES))

all.MSig <- read.gmt(file.path(data.dir,"ext","MSigDB.gmt.txt"))
all.MSig$ont <- toupper(all.MSig$ont)
all.MSig <- data.table(all.MSig)
m <- match(all.MSig$gene,row.names(df_tests_diff))
all.MSig$logRatio <- df_tests_diff$logRatio[m]
all.MSig$ttest.adj <- df_tests_diff$ttest.adj[m]

sig.GSEA.pos <- as.vector(GSEA_res[NES>0]$Description)
sig.GSEA.neg <- as.vector(GSEA_res[NES<0]$Description)

pos.MSig <- all.MSig[ont%in%sig.GSEA.pos]
pos.MSig.median <- pos.MSig[,list(med.z=mean(logRatio,na.rm=T),
							nb.gene.sig=length(gene[!is.na(logRatio) & logRatio>0 & ttest.adj>=2]),
								nb.gene=length(gene[!is.na(logRatio)])
								),by=ont]
pos.MSig.median$pc.sig <- pos.MSig.median$nb.gene.sig/pos.MSig.median$nb.gene
neg.MSig <- all.MSig[ont%in%sig.GSEA.neg]
neg.MSig.median <- neg.MSig[,list(med.z=mean(logRatio,na.rm=T),
							nb.gene.sig=length(gene[!is.na(logRatio) & logRatio<0 & ttest.adj>=2]),
								nb.gene=length(gene[!is.na(logRatio)])
								),by=ont]
neg.MSig.median$pc.sig <- neg.MSig.median$nb.gene.sig/neg.MSig.median$nb.gene

neg.group <- read.table(file.path(data.dir,"precomp","GSEA_neg_table_group.tab"),sep="\t")
pos.group <- read.table(file.path(data.dir,"precomp","GSEA_pos_table_group.tab.txt"),sep="\t")

m <- match(neg.MSig.median$ont,neg.group$V1)
neg.MSig.median$group_functions <- neg.group$V2[m]
neg.MSig.median <- na.omit(neg.MSig.median)

m <- match(pos.MSig.median$ont,pos.group$V1)
pos.MSig.median$group_functions <- pos.group$V2[m]
pos.MSig.median <- na.omit(pos.MSig.median)

MSig.median <- rbind(pos.MSig.median,neg.MSig.median)
MSig.median <- merge(MSig.median,GSEA_res,by.x="ont",by.y="Description")

write.table(MSig.median,file.path(res.dir,"TableSupp1C.tab"),sep="\t",quote=F,row.names=F)

group_functions <- levels(MSig.median$group_functions)
colour.functions <- c("cell_cycle"=colours()[34],"chromosomal_pos"=colours()[292],cinsarc="black",cancer="black",DNA_repair=colours()[367],DNA_replication=colours()[456],E2F_targets=colours()[33],MECOM_targets="grey",metabolism=colours()[31],mitochondria_activity=colours()[142],muscle=colours()[121],other="grey",RB1_targets=colours()[32],smooth_muscle=colours()[124],SRF_targets=colours()[131],collagen_organisation=colours()[57],EMT=colours()[11],ER_Golgi_related=colours()[84],ETS_family_targets=colours()[420],inflammatory_response=colours()[41],"mir-145_targets"=colours()[617],multicancer_invasiveness=colours()[367],MYC_targets=colours()[135],TGFB_signaling=colours()[8])

g <- ggplot(MSig.median,aes(x=med.z,y=abs(NES))) + 
	geom_point(aes(colour=group_functions,fill=group_functions,size=pc.sig),alpha=0.4) + 
	scale_colour_manual(values=colour.functions) + scale_fill_manual(values=colour.functions) + 
	theme(panel.background = element_rect(fill = "white", colour = "grey50"))  + 
    scale_size(range = c(2, 10)) + xlim(-4.5,4.5) + scale_shape_manual(values=c(functionnal=21,TFBS=24,miRNA=23,experimental=22)) +
    geom_text_repel(aes(label=ont),size=1)

pdf(file.path(res.dir,"Figure1D_GSEA_functional_annots.pdf"), paper="special",width=8,height=5)
print(g)
#export.plot(file.path(res.dir,"Figure1D_GSEA_functional_annots"),width=8,heigh=5)
dev.off()
message("Figure 1D OK")

####### iCistarget
message("iCistarget")

up.cistarget <- import.icistarget(dir=file.path(data.dir,"precomp","icistarget_hLMS_up"),group="up")
down.cistarget <- import.icistarget(dir=file.path(data.dir,"precomp","icistarget_hLMS_down"),group="down")

cistarget <- rbind(up.cistarget,down.cistarget)
cistarget[,fileName:=NULL]

cistarget$NES <- as.numeric(as.vector(cistarget$NES))
cistarget$weight <- (cistarget$nb.genes*cistarget$pc.ranked)

cistarget$weight.scaled <- (cistarget$weight/abs(max(cistarget$weight)))*5


cistarget.pwm <- cistarget[type=="PWMs"]
cistarget.pwm$weight <- as.numeric(as.vector(cistarget.pwm$weight))
cistarget.pwm$NES <- as.numeric(as.vector(cistarget.pwm$NES))
cistarget.pwm$nb.genes <- as.numeric(as.vector(cistarget.pwm$nb.genes))

cistarget.pwm <- cistarget.pwm[,list(weight=weight[which.max(abs(NES))],NES=NES[which.max(abs(NES))],nb.genes=nb.genes[which.max(abs(NES))],nb.dat=length(NES)),by=c("group","PWMs.cluster")]
cistarget.pwm <- cistarget.pwm[nb.dat>1]
cistarget.pwm <- cistarget.pwm[order(PWMs.cluster)][order(group)]
PWMs <- c("ZBTB42","XBP1","CREB3L1","CREB3L1","GC_rich_motif","SIP4","GC_rich_motif","unknown","CEBPB","ETS_family","SRF","plant","MEF2_family")

cistarget.pwm$TF <- PWMs
cistarget.pwm <- cistarget.pwm[,list(weight=weight[which.max(abs(NES))],NES=NES[which.max(abs(NES))],nb.genes=nb.genes[which.max(abs(NES))]),by=c("group","TF")]
cistarget.pwm <- cistarget.pwm[!TF%in%c("plant","fly","GC_rich_motif","unknown")]
cistarget.pwm$type <- "PWMs"
cistarget.pwm$tissue_celltype <- NA


cistarget.histones <- cistarget[type=="Histone modifications"]
cistarget.histones$TF <- sub(" ChIP-seq","",unlist(sapply(strsplit(cistarget.histones$Description," in "),"[[",1)))
cistarget.histones$tissue_celltype <- unlist(sapply(strsplit(as.vector(cistarget.histones$Description)," in "),"[[",2))
cistarget.histones$tissue_celltype <- sub(" $","",cistarget.histones$tissue_celltype)
#group          TF     weight      NES nb.genes type tissue_celltype
cistarget.histones <- cistarget.histones[,list(weight=weight[which.max(abs(NES))],NES=NES[which.max(abs(NES))],nb.genes=nb.genes[which.max(abs(NES))]),by=c("group","TF","type", "tissue_celltype")]
cistarget.histones <- cistarget.histones[,list(group,TF,weight,NES,nb.genes,type,tissue_celltype)]

cistarget.chip <- cistarget[type=="TF binding sites"]

Description <- gsub("ChIP-seq.*on human","",cistarget.chip$Description)
Description <- gsub("ECC-1. Note- This experiment previously referred to its biosample as ECC-1, however it has been found that all currently available ECC-1 are actually ","",Description)
Description <- gsub("ChIP-seq in","",Description)
Description <- gsub(" produced by the Snyder lab","",Description)
Description <- gsub("[ ]+"," ",Description)
Description <- strsplit(Description," ")
tissue_celltype <- sapply(Description,"[[",2)
TF <- c(sapply(Description[1:2],"[[",3),sapply(Description[-c(1:2)],"[[",1))
cistarget.chip$TF <- TF 
cistarget.chip$tissue_celltype <- tissue_celltype
cistarget.chip <- cistarget.chip[,list(weight=weight[which.max(abs(NES))],NES=NES[which.max(abs(NES))],nb.genes=nb.genes[which.max(abs(NES))]),by=c("group","TF","type", "tissue_celltype")]
cistarget.chip <- cistarget.chip[,list(group,TF,weight,NES,nb.genes,type,tissue_celltype)]

cistarget <- rbind(cistarget.chip,cistarget.histones,cistarget.pwm)
cistarget <- cistarget[abs(weight)>=3]
cistarget.group.feature <- rep(NA,length(cistarget$tissue_celltype))
cistarget.group.feature[grep("Smooth Muscle",cistarget$tissue_celltype)] <- "smooth_muscle"
cistarget.group.feature[grep("Aorta",cistarget$tissue_celltype)] <- "aorta"
cistarget.group.feature[grep("Ovary",cistarget$tissue_celltype)] <- "ovary"
cistarget.group.feature[grepl("Heart",cistarget$tissue_celltype) | grepl("Ventricle",cistarget$tissue_celltype)] <- "(fetal)_heart"
cistarget.group.feature[grep("IMR90",cistarget$tissue_celltype)] <- "IMR90"
cistarget.group.feature[cistarget$TF%in%c("CREB3L1","XBP1")] <- "UPR"
cistarget.group.feature[cistarget$TF%in%c("SRF","MEF2_family","MEF2A")] <- "muscle"
cistarget.group.feature[cistarget$TF%in%c("ZBTB42")] <- "skeletal_muscle"
cistarget.group.feature[grep("Derived Mesenchymal Stem Cell",cistarget$tissue_celltype)] <- "Derived MSC"
cistarget.group.feature[is.na(cistarget.group.feature) & cistarget$type!="Histone modifications"] <- cistarget$tissue_celltype[is.na(cistarget.group.feature) & cistarget$type!="Histone modifications"]
cistarget.group.feature[is.na(cistarget.group.feature)] <- "other"

cistarget$group.feature <- cistarget.group.feature

cistarget.up <- cistarget[group=="up"]
cistarget.up.hist <- cistarget.up[type=="Histone modifications"]
cistarget.up.hist$y <- as.numeric(factor(cistarget.up.hist$group.feature,levels=c("smooth_muscle","aorta","ovary","IMR90","(fetal)_heart","other")))
cistarget.up.TF <- cistarget.up[type!="Histone modifications"]
cistarget.up.TF$y <- as.numeric(factor(cistarget.up.TF$group.feature,levels=c("muscle","LoVo","IMR90","MCF-7","K562"))) + max(cistarget.up.hist$y) + 1

cistarget.dn <- cistarget[group=="down"]
cistarget.dn$NES <- -cistarget.dn$NES
cistarget.dn.hist <- cistarget.dn[type=="Histone modifications"]
cistarget.dn.hist$y <- as.numeric(factor(cistarget.dn.hist$group.feature,levels=c("Derived MSC","IMR90","other")))

cistarget.dn.TF <- cistarget.dn[type!="Histone modifications"]
cistarget.dn.TF$y <- as.numeric(factor(cistarget.dn.TF$group.feature,levels=c("UPR","skeletal_muscle","PFSK-1","SK-N-SH","K562","GM15510","other"))) + max(cistarget.up.hist$y) + 1

cistarget.plot <- rbind(cistarget.up.hist,cistarget.up.TF,cistarget.dn.hist,cistarget.dn.TF)
cistarget.plot$weight <- abs(cistarget.plot$weight)

cistarget.plot$celltype <- c("neuronal_Projenitor_CL", "fetal_lung_fibroblast_CL", "aorte", "aorte", "colon", "fetal_heart", "fetal_heart", "fetal_heart", "heart", "ovary", "ovary", "rectal", "rectal", "heart", "stomach", "stomach", "stomach", "CML_CL", "CML_CL", "endometrial_tumor", "ES_CL", "neuroblastoma_metastasis_CL", "fetal_lung_fibroblast_CL", "breast_tumor_CL", "colon_metastasis_CL", "prediction", "prediction", "MS_CL", "MS_CL", "fetal_lung_fibroblast_CL", "fetal_lung_fibroblast_CL", "fetal_lung_fibroblast_CL", "MS_CL", "breast", "skin_fibroblast", "monocyte_CD14", "astrocyte", "neuroectodermal_CL", "neuroblastoma_metastasis_CL", "neuroblastoma_metastasis_CL", "neuroblastoma_metastasis_CL", "neuroblastoma_metastasis_CL", "neuroblastoma_metastasis_CL", "neuroblastoma_metastasis_CL", "CML_CL", "lymphoblastoid_CL", "prediction","prediction","prediction","prediction","prediction","prediction")
celltype.color <- c(neuronal_Projenitor_CL=colours()[464],fetal_lung_fibroblast_CL=colours()[142],aorte=colours()[72],colon=colours()[121],colon_metastasis_CL=colours()[130],rectal=colours()[124],ovary=colours()[110],fetal_heart=colours()[59],heart=colours()[134],stomach=colours()[617],CML_CL=colours()[585],endometrial_tumor=colours()[640],ES_CL=colours()[517],MS_CL=colours()[518],neuroblastoma_metastasis_CL=colours()[99],breast=colours()[39],breast_tumor_CL=colours()[41],skin_fibroblast=colours()[148],monocyte_CD14=colours()[588],astrocyte=colours()[84],neuroectodermal_CL=colours()[544],lymphoblastoid_CL=colours()[584],prediction="grey")

g <- ggplot(cistarget.plot,aes(x=NES,y=y)) + 
	geom_segment(aes(x=0,xend=NES,y=y,yend=y),linetype="dashed",size=0.3,colour="grey") + 
	geom_vline(xintercept=0) + 	
	geom_point(aes(size=weight,shape=type,colour=celltype),alpha=0.4) +
	#geom_point(aes(size=weight,shape=type,colour=TF),alpha=0.4,data=cistarget.plot[type=="Histone_modifications"]) +
	theme(panel.background = element_rect(fill = "white", colour = "grey50")) +
	scale_size(range = c(2, 10)) + scale_colour_manual(values=celltype.color) +
	scale_y_reverse() + xlim(-16,16)

pdf(file.path(res.dir,"Figure1E_cistarget_results.pdf"), paper="special",width=8,height=6)
print(g)
#export.plot(file.path(res.dir,"Figure1E_cistarget_results"),width=8,heigh=6)
dev.off()
message("Figure 1E OK")

