# This script aims at analizing miRNA expression from TCGA cohort 
# From TCGA miRNA normalized data, apply limma for DE estimation and
# Holm's multi-testing correction
# from miRComb package
# miRNA source file : https://tcga.xenahubs.net/download/TCGA.SARC.sampleMap/miRNA_HiSeq_gene.gz - log2(RPM+1)
# mRNA source file : https://tcga.xenahubs.net/download/TCGA.SARC.sampleMap/HiSeqV2.gz
# Source necessary functions and libraries
message("Start miRNA analysis from TCGA")

source("3.functions.R")

# TCGA miRNA from SARC project - log2(RPM+1)
load(file.path(data.dir,"RData","TCGA_expr_miRNA.RData")) # miRNA_TCGA

annots.TCGA <- merge.clinical[merge.clinical$exp=="TCGA",]
# Select and order patients present in miRNA and mRNA data
miRNA_TCGA <- miRNA_TCGA[,colnames(miRNA_TCGA)%in%row.names(annots.TCGA) & colnames(miRNA_TCGA)%in%colnames(TCGA.expr)]
annots.TCGA <- annots.TCGA[row.names(annots.TCGA)%in%colnames(miRNA_TCGA),]
m <- match(row.names(annots.TCGA),colnames(miRNA_TCGA))
miRNA_TCGA <- miRNA_TCGA[,m]
m <- match(row.names(annots.TCGA),colnames(TCGA.expr))
TCGA.expr <- TCGA.expr[,m]

all.equal(row.names(annots.TCGA),colnames(TCGA.expr))
all.equal(row.names(annots.TCGA),colnames(miRNA_TCGA))
all.equal(colnames(miRNA_TCGA),colnames(TCGA.expr))

# Compute some statistics 
message("Compute some statistics on miRNA expression")

miRNA_TCGA[is.na(miRNA_TCGA)] <- 0
stats.expression <- data.table(t(miRNA_TCGA),group=annots.TCGA$Cluster)
median.expression <- stats.expression[,lapply(.SD,median),by=group]
median.expression.miRNA <- as.data.frame(t(median.expression[,group:=NULL]))
colnames(median.expression.miRNA) <- c("median.hLMS","median.oLMS")
median.expression.miRNA$ID <- row.names(median.expression.miRNA)
median.expression.miRNA <- data.table(median.expression.miRNA)


message("run mirComb")
pheno.data <- data.frame(group=annots.TCGA$Cluster,DvH=ifelse(annots.TCGA$Cluster=="hLMS",1,0),row.names=row.names(annots.TCGA))
data.obj<-new("corObject",dat.miRNA=as.matrix(miRNA_TCGA),dat.mRNA=as.matrix(TCGA.expr),
                           pheno.miRNA=pheno.data,pheno.mRNA=pheno.data)
data.obj<-addDiffexp(data.obj,"miRNA",classes="DvH",method.dif="limma",method.adj="holm")
data.obj<-addDiffexp(data.obj,"mRNA",classes="DvH",method.dif="limma",method.adj="holm")

# extract DE objects

miRNA.diffexp <- data.obj@diffexp.miRNA
miRNA.diffexp$ID <- row.names(miRNA.diffexp)
miRNA.diffexp <- data.table(miRNA.diffexp)
miRNA.diffexp <- merge(miRNA.diffexp,median.expression.miRNA,by="ID")

# Figure 2A : TCGA PCA
# Select miRNAs with a median of at least 1 normalized count in at least one group
message("Compute PCA on all expressed miRNA expression")
miRNA_TCGA <- miRNA_TCGA[row.names(miRNA_TCGA)%in%miRNA.diffexp[median.hLMS>1 | median.oLMS>1]$ID,]
pca.LMS <- prcomp(t(miRNA_TCGA))

g <- ggbiplot(pca.LMS,choices=c(1,2), obs.scale = 1, var.scale = 1, 
                      groups =annots.TCGA$Cluster[na.omit(match(colnames(miRNA_TCGA),row.names(annots.TCGA)))], 
                      ellipse = T, var.axes=F,
                      circle = F,varname.size = 0) +
  theme(legend.direction = 'horizontal', legend.position = 'top',aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black")) +
  scale_colour_manual(values=c("hLMS"=colours()[613],"oLMS"=colours()[132])) 
pdf(file.path(res.dir,"Figure2A_PCA_TCGA_miRNA.pdf"),width=4,heigh=4)
print(g)
#export.plot(file.path(res.dir,"Figure2A_PCA_TCGA_miRNA"),width=4,heigh=4)
dev.off()
message("Figure 2A TCGA OK ")

# Figure 2C : TCGA heatmap

message("Heatmap on differentially expressed miRNA expression")
miRNA_TCGA_diff <- miRNA_TCGA[row.names(miRNA_TCGA)%in%miRNA.diffexp[adj.pval<0.01 & abs(logratio)>1]$ID,]
mat <- miRNA_TCGA_diff
mat[is.na(mat)] <- 0
res.pca <- PCA(t(mat),scale.unit=FALSE,ncp=2, graph = FALSE)
res.hcpc <- HCPC(res.pca, nb.clust=4, conso=0,graph=F,order=F)

mat.temp <- mat[,match(res.hcpc$call$t$tree$labels,colnames(mat))]

res.pca <- PCA(mat.temp,scale.unit=FALSE,ncp=10, graph = FALSE)
res.hcpc.row <- HCPC(res.pca, nb.clust=4, conso=0,graph=F,order=F)

mat.temp <- mat.temp[match(res.hcpc.row$call$t$tree$labels,row.names(mat.temp)),]
mat.temp <- t(apply(mat.temp,1,function(x){(x - mean(x))/sd(x)}))

pdf(file.path(res.dir,"Figure2C_heatmap_diff_miRNA_cr_TCGA.pdf"),width=5,heigh=3)
pheatmap(mat.temp,annotation_col=data.frame(Cluster=merge.clinical$Cluster[merge.clinical$Cluster!="gLMS"],row.names=row.names(merge.clinical)[merge.clinical$Cluster!="gLMS"]),
         color=colorRampPalette(c(colours()[131],colours()[131],"white","red","red"))(100),breaks=seq(-5,5,0.1),#cluster_rows=res.hcpc.row$call$t$tree,
         annotation_colors=list(Cluster=c("hLMS"=colours()[613],"oLMS"=colours()[132])),show_rownames=F, cluster_cols=res.hcpc$call$t$tree,
         border_color=NA, clustering_distance_row="correlation",
         show_colnames=F)
#export.plot(file.path(res.dir,"Figure2C_heatmap_diff_miRNA_cr_TCGA"),width=5,heigh=3)
dev.off()
message("Figure 2C TCGA OK ")








