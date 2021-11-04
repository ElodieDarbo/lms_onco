message("Start GTEX analysis with hLMS top 100 highest expressed genes")


message("loading GTEX TCGA data")
load(file.path(data.dir,"RData","GTEX_objetcs.RData")) # GTEX_TCGA,samples.all,batch,tissue
message("Compute hLMS median expression")
hLMS <- apply(GTEX_TCGA[,samples.all$Cluster=="hLMS"],1,median)
hLMS <- sort(hLMS,decreasing=T)
message("Filter GTEX object ")
select.GTEX <- GTEX_TCGA[row.names(GTEX_TCGA)%in%names(hLMS)[1:100],]


set.seed(100)
message("Computing tSNE clustering")
tsne.group <- Rtsne(t(select.GTEX), dims = 2, perplexity=150, verbose=TRUE, max_iter = 1500)
tsne.group <- as.data.frame(tsne.group$Y)
colnames(tsne.group) <- c("tSNE1","tSNE2")
tsne.group$samples <-  samples.all$Cluster[match(colnames(select.GTEX),row.names(samples.all))]
tsne.group$LMS <- ifelse(tsne.group$samples%in%c("hLMS","oLMS"),as.vector(tsne.group$samples),"other")
g <- ggplot(aes(x=tSNE1, y=tSNE2, group=samples), data=tsne.group) + 
  geom_point(aes(x=tSNE1, y=tSNE2, fill=samples,color=samples),shape=21,alpha=1,stroke=0.1,size=ifelse(tsne.group$LMS=="other",1,1.8)) +   
  stat_ellipse(aes(color=samples),data=tsne.group[tsne.group$LMS!="other",]) +
  theme(panel.background = element_rect(fill = "white", colour = "black") , legend.position="right",aspect.ratio=1) +
  geom_hline(yintercept=0,linetype="dashed",size=0.3) +
  geom_vline(xintercept=0,linetype="dashed",size=0.3) +
  scale_color_manual(values=tissue) + scale_fill_manual(values=tissue) 

pdf(file.path(res.dir,"FS2A_tSNE_top100_hLMS_GTEX.pdf"),width=8,heigh=12)
plot(g)
dev.off()
#export.plot(file.path(res.dir,"FS2A_tSNE_top100_hLMS_GTEX"),width=12,heigh=8)
message("FS2A OK")

message("Computing median expression per tissue")
median.tissue <- do.call("cbind",lapply(unique(samples.all$Cluster),function(x,all.TCGA,samples.all){
  hLMS <- apply(all.TCGA[,samples.all$Cluster==x],1,median)
},select.GTEX,samples.all))
colnames(median.tissue) <- unique(samples.all$Cluster)

pdf(file.path(res.dir,"FS2B_heatmap_top100_median_tissue.pdf"),width=14,heigh=7)
pheatmap(median.tissue,border_color = NA,
         cluster_rows = T,color = colorRampPalette(c("white","white","white",colours()[404],colours()[404],"red","brown","brown","black"))(100),
         clustering_distance_cols = "correlation",clustering_distance_rows = "correlation",clustering_method = "average")
#export.plot(file.path(res.dir,"FS2B_heatmap_top100_median_tissue"),heigh=14, width=7)
dev.off()
message("FS2B OK")

