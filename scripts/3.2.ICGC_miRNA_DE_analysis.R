# This script aims at analysis miRNA data from both ICGC
# ICGC: from raw counts data, using edgeR normalsation and differential expression (DE)
# fitting a generalized linear model. 
message("Start miRNA analysis from ICGC")

annots.ICGC <- merge.clinical[merge.clinical$exp=="ICGC",]


message("Computing edgeR differential expression")
# load miRNA raw counts
load(file.path(data.dir,"RData","ICGC_miRNA_raw_counts.RData")) # raw.ICGC.miRNA

group <- factor(annots.ICGC$Cluster[match(colnames(raw.ICGC.miRNA),row.names(annots.ICGC))],levels=c("oLMS","hLMS"))
design <- model.matrix(~group)

edger.list <- DGEList(counts=raw.ICGC.miRNA, group=group)

keep <- filterByExpr(edger.list,min.count = 5, min.total.count = 10)
counts.filt <-  raw.ICGC.miRNA[keep,]
edger.filt <- DGEList(counts=counts.filt, group=group)

edger.norm <- calcNormFactors(edger.filt)
edger.norm <- estimateDisp(edger.norm,design)

miRNA_ICGC_norm <- log2(cpm(edger.norm)+1)

fit <- glmQLFit(edger.norm)
qlf.2vs1 <- glmQLFTest(fit, coef=2)

edger.GLM <- as.data.frame(topTags(qlf.2vs1,n = Inf))
edger.GLM$ID <- row.names(edger.GLM)
edger.GLM <- data.table(edger.GLM)
edger.GLM$adj.pval <- p.adjust(edger.GLM$PValue,method = "holm")

message("Compute some statistics on miRNA expression")
median.hLMS <- apply(miRNA_ICGC_norm[,group=="hLMS"],1,median)
median.oLMS <- apply(miRNA_ICGC_norm[,group=="oLMS"],1,median)
median.ICGC <- data.frame(median.hLMS,median.oLMS)
median.ICGC$ID <- row.names(median.ICGC)
median.ICGC <- data.table(median.ICGC)

miRNA.ICGC <- merge(edger.GLM,median.ICGC,by="ID")
#save(miRNA_ICGC_norm,miRNA.ICGC,file=file.path(res.dir,"ICGC_miRNA_norm_counts_DE.RData"))

#write.table(miRNA.ICGC,file.path(res.dir,"ICGC_miRNA_DE.tab"),sep="\t",row.names=F,quote=F)


# Figure 2A : ICGC PCA
message("Compute PCA on all expressed miRNA expression")
pca.ICGC <- prcomp(t(miRNA_ICGC_norm))

g <- ggbiplot(pca.ICGC,choices=c(1,2), obs.scale = 1, var.scale = 1, 
                   groups = annots.ICGC$Cluster[na.omit(match(colnames(miRNA_ICGC_norm),row.names(annots.ICGC)))], 
                   ellipse = T, var.axes=F,
                   circle = F,varname.size = 0) +
  theme(legend.direction = 'horizontal', legend.position = 'top',aspect.ratio=1,panel.background = element_rect(fill = "white", colour = "black")) +
  scale_colour_manual(values=c(oLMS=colours()[132],hLMS=colours()[613]))
pdf(file.path(res.dir,"Figure2A_PCA_ICGC_miRNA.pdf"),width=4,heigh=4)
print(g)
#export.plot(file.path(res.dir,"Figure2A_PCA_ICGC_miRNA"),width=8,heigh=4)
dev.off()
message("Figure 2A ICGC OK ")

# Figure 2C : ICGC heatmap
message("Heatmap on differentially expressed miRNA expression")
miRNA_ICGC_diff <- miRNA_ICGC_norm[row.names(miRNA_ICGC_norm)%in%miRNA.ICGC[adj.pval<0.01 & abs(logFC)>1]$ID,]
mat <- miRNA_ICGC_diff
mat[is.na(mat)] <- 0
res.pca <- PCA(t(mat),scale.unit=T,ncp=2, graph = FALSE)
res.hcpc <- HCPC(res.pca, nb.clust=4, conso=1,graph=F,order=F)

mat.temp <- mat[,match(res.hcpc$call$t$tree$labels,colnames(mat))]
mat.temp <- t(apply(mat.temp,1,function(x){(x - mean(x))/sd(x)}))

pdf(file.path(res.dir,"Figure2C_heatmap_diff_miRNA_centered_reduced_ICGC.pdf"),width=5,heigh=3)

pheatmap(mat.temp,annotation_col=data.frame(Cluster=merge.clinical$Cluster[merge.clinical$Cluster!="gLMS"],row.names=row.names(merge.clinical)[merge.clinical$Cluster!="gLMS"]),
         color=colorRampPalette(c(colours()[131],colours()[131],"white","red","red"))(100),breaks=seq(-4,4,0.08),#cluster_rows=res.hcpc.row$call$t$tree,
         annotation_colors=list(Cluster=c("hLMS"=colours()[613],"oLMS"=colours()[132])),show_rownames=F, cluster_cols=res.hcpc$call$t$tree,
         border_color=NA, clustering_distance_row="correlation",
         show_colnames=F)
#export.plot(file.path(res.dir,"Figure2C_heatmap_diff_miRNA_centered_reduced_ICGC"),width=5,heigh=3)
dev.off()
message("Figure 2C ICGC OK ")



