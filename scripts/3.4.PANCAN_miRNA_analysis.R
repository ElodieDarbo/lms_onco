message("Start TCGA PANCAN miRNA expression analysis")


##########################################################
##########################################################
# TCGA PANCAN clinical and expression data
load(file.path(data.dir,"RData","clinical_miRNA_PANCAN.RData")) # clinical.miRNA.PANCAN
clinical <- clinical.miRNA.PANCAN

# Import colors for cancer types
X_primary_disease <- import.PANCAN.colors()
# TCGA PANCAN batch corrected miRNA normalized expression 
miRNA.TCGA <- read.table(file.path(data.dir,"ext","batch_corrected_miRNA_PANCAN.tab"),sep="\t",head=T,row=1)
m <- match(colnames(miRNA.TCGA),clinical$sampleID)
miRNA.TCGA <- miRNA.TCGA[,!is.na(m)]
clinical <- clinical[na.omit(m),]
##########################################################
##########################################################

##########################################################
##########################################################
message("Heatmap on differentially expressed miRNA in both cohorts")
# differential miRNA
diff.miRNA <-  miRNA_ICGC_TCGA[Significance=="TCGA+ICGC"]$ID
# Add mir-1 because hsa-miR-1-3p is not in PANCAN
miRNA.TCGA.diff <- miRNA.TCGA[row.names(miRNA.TCGA)%in%c(diff.miRNA,"hsa-miR-1"),]
miRNA.TCGA.diff[is.na(miRNA.TCGA.diff)] <- 0
mat <- miRNA.TCGA.diff

res.pca <- PCA(t(mat),scale.unit=F,ncp=3, graph = FALSE)
res.hcpc <- HCPC(res.pca, nb.clust=11, conso=0,graph=F,order=F)

temp.cluster <- data.frame(sampleID=row.names(res.hcpc$data.clust),Cluster=res.hcpc$data.clust$clust)
temp.cluster$X_primary_disease <- clinical$X_primary_disease[match(temp.cluster$sampleID,clinical$sampleID)]

message("v test on cluster containing hLMS")
stats.hcpc <- res.hcpc$desc.var$quanti$`4`
#write.table(stats.hcpc,file.path(res.dir,"cluster_miRNA_PANCAN_discriminative_hLMS_HCPC.tab"),sep="\t",quote=F)
print(stats.hcpc)
mat.temp <- mat[,match(res.hcpc$call$t$tree$labels,colnames(mat))]
mat.temp <- mat.temp - apply(mat.temp,1,median)

pdf(file.path(res.dir,"Figure2B_heatmap_diff_miRNA_all_TCGA.pdf"),width=10,heigh=5)

clust <- pheatmap(mat.temp,show_rownames=F,show_colnames=F,annotation_col=data.frame(X_primary_disease=clinical$X_primary_disease,row.names=clinical$sampleID),
                  cutree_cols=11, cluster_cols=res.hcpc$call$t$tree, border_color=NA,
                  annotation_colors=list(X_primary_disease=X_primary_disease),
                  color=colorRampPalette(c(colours()[132],colours()[121],"white","red","brown"))(100),breaks=seq(-11,11,0.22))
#export.plot(file.path(res.dir,"Figure2B_heatmap_diff_miRNA_all_TCGA"),width=10,heigh=5)
dev.off()
message("Figure2B heatmap OK")

cluster=factor(cutree(clust$tree_col, k = 11))
clinical$PANCAN.diff <- cluster

message("Barplot distribution cancer type after HCPC clustering")

melt.PANCAN <- as.data.frame.matrix(table(clinical$X_primary_disease,clinical$PANCAN.diff))
melt.PANCAN$type <- row.names(melt.PANCAN)
melt.PANCAN <- melt(melt.PANCAN)
melt.PANCAN$variable <- factor(as.vector(melt.PANCAN$variable),levels=c(1,9,3,11,10,7,2,6,8,5,4))
g <- ggplot(melt.PANCAN,aes(x=variable,y=value,fill=type)) + 
  geom_col(position="fill") + scale_fill_manual(values=X_primary_disease) +
  theme(legend.position="bottom",panel.background = element_rect(fill = "white", colour = "black") )
pdf(file.path(res.dir,"Figure2B_PANCAN_cluster_diff_histotype.pdf"),width=14,heigh=7)
print(g)
#export.plot(file.path(res.dir,"Figure2B_PANCAN_cluster_diff_histotype"),width=14,heigh=7)
dev.off()
message("Figure2B barplot OK")

message("Count number of samples per cancer type per clusters")
counting.samples <- as.data.frame.matrix(table(clinical$X_primary_disease,clinical$PANCAN.diff))
counting.samples$total <- rowSums(counting.samples)
print(counting.samples)
##########################################################
##########################################################
# PANCAN DIO3 genes

message("Heatmap on DIO-DLK1 miRNA cluster")
hairpins.pancan <- tolower(sub("-[35]p","",row.names(miRNA.TCGA)))
miRNA.TCGA.DIO3 <- miRNA.TCGA[hairpins.pancan%in%DIO3,]
miRNA.TCGA.DIO3[is.na(miRNA.TCGA.DIO3)] <- 0

# Figure supp 3
# heatmap
mat <- miRNA.TCGA.DIO3

res.pca <- PCA(t(mat),scale.unit=F,ncp=3, graph = FALSE)
res.hcpc <- HCPC(res.pca, nb.clust=11, conso=0,graph=F,order=F)

mat.temp <- mat[,match(res.hcpc$call$t$tree$labels,colnames(mat))]
#mat.temp <- mat.temp - apply(mat.temp,1,median)

pdf(file.path(res.dir,"FS3C_heatmap_DIO3_miRNA_all_TCGA.pdf"),width=10,heigh=5)

clust <- pheatmap(mat.temp,show_rownames=F,show_colnames=F,annotation_col=data.frame(X_primary_disease=clinical$X_primary_disease,row.names=clinical$sampleID),
                  cutree_cols=12, cluster_cols=res.hcpc$call$t$tree, border_color=NA,
                  annotation_colors=list(X_primary_disease=X_primary_disease),
                  color=colorRampPalette(c("white","white",colours()[420],"red","red","brown"))(100))

#export.plot(file.path(res.dir,"FS3C_heatmap_DIO3_miRNA_all_TCGA"),width=10,heigh=5)
dev.off()
message("FigureS3C heatmap OK")

# Use clustering from HCPC
message("Barplot distribution cancer type after HCPC clustering")

cluster=factor(cutree(clust$tree_col, k = 12))
clinical$PANCAN.DIO3 <- cluster

# Barplot number of samples per cancer type per clusters
melt.PANCAN.DIO3 <- as.data.frame.matrix(table(clinical$X_primary_disease,clinical$PANCAN.DIO3))
melt.PANCAN.DIO3$type <- row.names(melt.PANCAN.DIO3)
melt.PANCAN.DIO3 <- melt(melt.PANCAN.DIO3)
melt.PANCAN.DIO3$variable <- factor(as.vector(melt.PANCAN.DIO3$variable),levels=c(10,9,7,5,8,1,3,11,12,4,6,2))
g <- ggplot(melt.PANCAN.DIO3,aes(x=variable,y=value,fill=type)) + 
  geom_col(position="fill") + scale_fill_manual(values=X_primary_disease) +
  theme(legend.position="bottom",panel.background = element_rect(fill = "white", colour = "black") )
pdf(file.path(res.dir,"FS3C_PANCAN_cluster_DIO3_histotype.pdf"),width=14,heigh=7)
print(g)
#export.plot(file.path(res.dir,"FS3C_PANCAN_cluster_DIO3_histotype"),width=14,heigh=7)
dev.off()
message("FigureS3C barplot OK")


message("Computing median expression of DIO3-DLK1 miRNA per patient")
# DIO3-DLK1 median expression per sample
mean.expr.sample <- apply(mat.temp,2,mean)
test.mean.expr <- data.table(mean.expr.sample,cluster,LMS=clinical$X_primary_disease)
test.mean.expr$LMS <- ifelse(test.mean.expr$LMS%in%c("hLMS","oLMS"),test.mean.expr$LMS,"other")
test.mean.expr$cluster <- factor(as.vector(test.mean.expr$cluster),levels=c(10,9,7,5,8,1,3,11,12,4,6,2))
test.mean.expr <- test.mean.expr[order(cluster)]

g <- ggplot(test.mean.expr,aes(y=mean.expr.sample,x=1:nrow(test.mean.expr))) + 
  geom_point(colour=colours()[231],size=0.3) +
  geom_point(aes(color=LMS),size=0.5) +
  scale_colour_manual(values=c("hLMS"=colours()[613],"oLMS"=colours()[132],other="NA")) +
  #geom_smooth(colour="black",size=0.3,fill=NA,formula = y ~ poly(x,9) ,method="lm") +
  theme(panel.background = element_rect(fill = "white", colour = "black"))
pdf(file.path(res.dir,"FigureS3C_points_cluster_DIO3_DLK1_all_TCGA.pdf"),width=14,heigh=7)
print(g)
#export.plot(file.path(res.dir,"FigureS3C_points_cluster_DIO3_DLK1_all_TCGA"),width=8,heigh=3)
dev.off()
message("FigureS3C points OK")

message("Count number of samples per cancer type per clusters")
counting.samples.DIO3 <- as.data.frame.matrix(table(clinical$X_primary_disease,clinical$PANCAN.DIO3))
counting.samples.DIO3$total <- rowSums(counting.samples.DIO3)
print(counting.samples.DIO3)







