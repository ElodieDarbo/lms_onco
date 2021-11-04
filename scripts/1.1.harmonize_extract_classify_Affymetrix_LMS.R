# This script aims at identify genes with enough replicability between 
# Agilent and Affymetrix platform through 87 patients analyzed on both
# platforms and harmonized the data to be able to combine them. 
# All datasets are then merged and harmonized. LMSs grouped in a single cluster
# are compared to others. Definition of hLMS and oLMS in Affymetrix cohort.

message("Start affymetrix cohort analysis")

# Source necessary functions and libraries
source("1.functions.R")

##########################################################################################
# PREPARE DATA 87 samples

message("Load micro-array data, pre-treated gcRMA + quantile normalized (1st round)")

load(file.path(data.dir,"RData","datasetlist.RData"))
# Agilent
agilent.87 <- datasetlist$sarcompAgi.87
colnames(agilent.87) <- sub(".sarcComp","",colnames(agilent.87))
# Sample IDs
samples <- colnames(agilent.87)

agilent.87 <- agilent.87[,samples]



# Load pre-processed data from RData objects
# Affymetrix
affy.389 <- datasetlist$sarcompAffy.389
colnames(affy.389) <- sub(".affy.sarcComp","",colnames(affy.389))
affy.87 <- affy.389[,samples]


message("Selecting genes with consistent gene profiles in 87 patients analysed on Affymetrix and Agilent platforms")
# Probes
# Affymetrix
affy.probes <- format.affy.probes(dir=file.path(data.dir,"ext"))

# Agilent
agilent.probes <- read.table(file.path(data.dir,"ext","Agilent.txt"),sep="\t",head=TRUE)
agilent.probes <- agilent.probes[agilent.probes[,2]!="",] # 30936 probes # 19595 gene symbols

## select genes represented in both datasets
all.genes <- as.vector(unique(agilent.probes$GeneSymbol[agilent.probes$GeneSymbol%in%as.vector(affy.probes$gene.symbol)]))
list.common <- list(affymetrix=unique(as.vector(affy.probes$gene.symbol)),agilent=unique(as.vector(agilent.probes$GeneSymbol)))
venn.all <- Venn(list.common)
pdf(file.path(res.dir,"FigureS1B.pdf"), paper="special",width=4,height=4)
plot(venn.all)
dev.off()
#export.plot(file.path(res.dir,"FigureS1B"),width=4,heigh=4)
message("Plot Figure S1B OK")

##########################################################################################
# order genes
m <- match(row.names(agilent.87),row.names(affy.87))
agilent.87 <- agilent.87[!is.na(m),]
affy.87 <- affy.87[na.omit(m),]

# verify samples are ordered the same way
all.equal(colnames(affy.87),colnames(agilent.87))

# This computation lasts around 1h30, you can either load the pre-computed object 

if (run.preProcessed){
	load(file.path(data.dir,"precomp","correlation_all_genes_87_samples.RData"))
}else{
	stats.all.genes <- stat.probes(affy.87,agilent.87)
	save(stats.all.genes,file=file.path(res.dir,"correlation_all_genes_87_samples.RData"))
}
# 

message("Select consistent genes")
# Select genes with better correlation between them than with any other gene in any platform
keep.genes <- select.genes.by.cor(stats.all.genes)
# save the genes to be use during harmonization
# save(keep.genes,file=file.path(res.dir,"selected_genes_correlating_platforms.RData"))

# verify genes are ordered the same way
all.equal(row.names(affy.87),row.names(agilent.87))

# create a merged data.frame
affy.agilent <- cbind(affy.87,agilent.87)
affy.agilent <- affy.agilent[row.names(affy.agilent)%in%keep.genes,]
colnames(affy.agilent) <- c(paste0("samp",1:87,".affy"),paste0("samp",1:87,".agi"))

message("Quantile normalization on merged datasets (2n round)")
## quartiles normalisation
all.norm <- normalize.quantiles(as.matrix(affy.agilent))
row.names(all.norm) <- row.names(affy.agilent)
colnames(all.norm) <- colnames(affy.agilent)

message("Harmonize by the medians")
## harmonised by medians
harmonized <- apply(all.norm,1,harmonize,n=87)
harmonized <- t(harmonized)
row.names(harmonized) <- row.names(all.norm)
colnames(harmonized) <- colnames(all.norm)

# FigureS1C Plot an example of normalization effects
# Line plots
pos1 <- 1
pos2 <- 50
g.raw <- plot.couple.gene(affy.agilent,pos1,pos2,"norm.by.tech")
g.norm1 <- plot.couple.gene(all.norm,pos1,pos2,"all.norm")
g.norm2 <- plot.couple.gene(harmonized,pos1,pos2,"hamonized")

pdf(file.path(res.dir,"FS1_example_transformation_87.pdf"), paper="special",width=8,height=7)
multiplot(g.raw,g.norm1,g.norm2,cols=1)
dev.off()
#export.plot(file.path(res.dir,"FS1_example_transformation_87"),width=8,heigh=7)

# Heatmaps plots
colors.heat <- colorRampPalette(c("yellow","orange","red"))(100)

# select genes varying in harmonized data to visualize all steps
var.genes <- apply(harmonized,1,var)
mat <- affy.agilent[var.genes>3,]
pdf(file.path(res.dir,"FS1_raw_heatmap_87_samples.pdf"), paper="special",width=6,height=7)
pheatmap(mat,show_colnames=F,show_rownames=F,color=colors.heat)
#export.plot(file.path(res.dir,"FS1_raw_heatmap_87_samples"),width=6,heigh=7)
dev.off()

mat <- all.norm[var.genes>3,]
pdf(file.path(res.dir,"FS1_Qnorm_heatmap_87_samples.pdf"), paper="special",width=6,height=7)
pheatmap(mat,show_colnames=F,show_rownames=F,color=colors.heat)
#export.plot(file.path(res.dir,"FS1_Qnorm_heatmap_87_samples"),width=6,heigh=7)
dev.off()

mat <- harmonized[var.genes>3,]
pdf(file.path(res.dir,"FS1_harmonized_heatmap_87_samples.pdf"), paper="special",width=6,height=7)
pheatmap(mat,show_colnames=F,show_rownames=F,color=colors.heat)
#export.plot(file.path(res.dir,"FS1_harmonized_heatmap_87_samples"),width=6,heigh=7)
dev.off()
message("Plots Figure S1C OK")


##########################################################################################
# HARMONIZE DATA all datasets
# prepare a merged datasets 
message("Harmonize all datasets")
raw.all <- lapply(c("sarcompAffy.389","sarcompAgi.87","lignees.39","LPS.50","SSX.58","GIST.60"),function(x,datasetlist,keep.genes){
	x <- datasetlist[[x]]
	m <- match(keep.genes,row.names(x))
	x <- x[na.omit(m),]
	return(x)
	},datasetlist,keep.genes)

raw.all <- do.call("cbind",raw.all)

#first round of normalisation
all.samples.norm <- as.data.frame(normalize.quantiles(as.matrix(raw.all)))
row.names(all.samples.norm) <- row.names(raw.all)
colnames(all.samples.norm) <- colnames(raw.all)

# Harmonization
harmonized.all.samples <- apply(all.samples.norm,1,harmonize.2,n=c(389,(87+39+50+58+60)))
harmonized.all.samples <- t(harmonized.all.samples)
row.names(harmonized.all.samples) <- row.names(all.samples.norm)
colnames(harmonized.all.samples) <- colnames(all.samples.norm)

# Remove cell line and duplicated samples
harmonized.all.samples <- harmonized.all.samples[,!grepl("MFH148LP83|MFH100P33",colnames(harmonized.all.samples)) & grepl("affy|GIST$|SSX_IB$|LPS$",colnames(harmonized.all.samples))]
colnames(harmonized.all.samples) <- sub(".affy.sarcComp|.GIST$|.SSX_IB$|.LPS$","",colnames(harmonized.all.samples) )

##########################################################################################
# Find co-expressed genes in the
# combined dataset
message("Find co-expressed genes in the combined dataset")

#load(file.path(res.dir,"harmonized_all_samples.RData"))

#save(harmonized.all.samples,file=file.path(res.dir,"harmonized_all_samples.RData"))
var.genes <- apply(harmonized.all.samples,1,var)
mat <- harmonized.all.samples[var.genes>2,]

mat <- (mat - apply(mat,1,median))/apply(mat,1,sd)

# Compute gene pairwise correlation and select highly correlated genes
dissimilarity <- cor(t(mat))
diag(dissimilarity) <- 0
dissimilarity[lower.tri(dissimilarity)] <- 0
dissimilarity <- as.data.frame(dissimilarity)
dissimilarity$int2 <- row.names(dissimilarity)
dissimilarity <- melt(dissimilarity)
colnames(dissimilarity) <- c("int2","int1","cor")
cor.test <- dissimilarity[dissimilarity$cor>=0.7,]

# produce a graph and search for communities
g.cor <- graph.data.frame(cor.test, directed=FALSE)
ebc.cor <- edge.betweenness.community(g.cor, directed=F)
coms.cor <- communities(ebc.cor)
new.modules <- coms.cor
new.modules <- new.modules[order(unlist(lapply(new.modules,length)),decreasing=T)]
new.modules <- new.modules[lapply(new.modules,length)>1]

# import pathways from MSig database
pathways <- import.pathways(data.dir)

# compute functional enrichment for modules with at least 5 genes

message("Annotate co-expressed gene modules")
i <- 1
modules.functions <- lapply(new.modules[unlist(lapply(new.modules,length))>=5],function(gene,pathways){
	m.tmp <- rbindlist(lapply(pathways,function(p,gene){
		egmt <- as.data.frame(enricher(gene,TERM2GENE=p))
		},gene))
	m.tmp <- m.tmp[order(p.adjust)]
	write.table(m.tmp,file.path(res.dir,paste("enrichment_cluster_multiplatform_",i,sep="")),sep="\t",row.names=F,quote=F)	
	i <<- i + 1
	return(m.tmp)
	},pathways=pathways[1:5])

#save(new.modules,modules.functions,file=file.path(res.dir,"modules_multiplatform_objects.RData"))

##########################################################################################
# Find co-expressed genes in the
# combined dataset

clinical <- import.clinical(dir=file.path(data.dir,"annotations"),names.samp=colnames(harmonized.all.samples))

# create a dataframe in order to annotate the heatmap with module enriched pathways
i <- 1
df.row <- rbindlist(lapply(new.modules,function(x){
	if (length(x)<5){
		res <- data.frame(SYMBOL=x,cluster=0)
	}
	else {
		res <- data.frame(SYMBOL=x,cluster=i)
		i <<- i + 1
	}
	res
	}))
df.row <- data.frame(modules=factor(df.row$cluster),row.names=df.row$SYMBOL)
df.row$modules <- as.factor(df.row$modules)
# This is manual annotation from interpretation of enrichment_cluster_multiplatform* files 
# created previously and kept in results folder
levels(df.row$modules) <- c("modules < 5 (89)","1. immune system activation (117)", "2. cell cycle (89)","3. leukocyte (32)",
	"4. smooth muscle related (21)","5. skeletal muscle related (18)", "6. viral response (12)", 
	"7. Adipogenesis (9)", "8. chr12q13-14 (8)", "9. interferon gamma response (7)", "10. NE (7)", "11. extracellular matrix (6)", "12. muscle contraction (6)",
	"13. chryq11-p11 (6)", "14. apical plasma membrane (5)","15. NE (5)")

# Format clinical data
df.cols <- clinical[,c("differentiation","mitotic_count_10_HPF","grade_primary_tumour","sex","Localisation","Type")]
df.cols$differentiation <- as.factor(df.cols$differentiation)
df.cols$grade_primary_tumour <- as.factor(df.cols$grade_primary_tumour)
df.cols$Localisation <- as.vector(df.cols$Localisation)
df.cols$sex <- as.vector(df.cols$sex)

# select and order genes according to genes present in co-expression modules
m <- match(row.names(df.row),row.names(harmonized.all.samples))
mat <- harmonized.all.samples[m,]

# FactoMineR clustering
res.pca <- PCA(t(mat),scale.unit=FALSE, ncp=14, graph = FALSE)
res.hcpc.harmonized <- HCPC(res.pca, nb.clust=5, conso=0, graph = FALSE)
df.cols$Cluster <- res.hcpc.harmonized$data.clust$clust[match(row.names(df.cols),row.names(res.hcpc.harmonized$data.clust))]

# ordering matrix column to be displayed with HCPC clustering results
mat.temp <- mat[,match(res.hcpc.harmonized$call$t$tree$labels,colnames(mat))]

# scaled data for displaying
harmonized.scaled <- (mat.temp - apply(mat.temp,1,median)) / apply(mat.temp,1,sd)


# Extraction on Leiomyosarcomas and prelimininary grouping hLMS/oLMS
LMS.select <- row.names(df.cols)[!is.na(df.cols$Type) & df.cols$Type=="leiomyosarcomes"]
df.LMS.extract <- df.cols[row.names(df.cols)%in%LMS.select,]
df.LMS.extract$Cluster <- ifelse(df.LMS.extract$Cluster==2,"hLMS","oLMS")

# Create a LMS dataset
affy.389 <- datasetlist$sarcompAffy.389
colnames(affy.389) <- sub(".affy.sarcComp","",colnames(affy.389))
puces.expr <- affy.389[,row.names(df.LMS.extract)]

#save(df.LMS.extract,puces.expr,file=file.path(res.dir,"LMS_extraction_objects.RData"))

##########################################################################################
# Compute differential expression
# 
message("Compute statistics hLMS/oLMS")
ttest.test <- apply(puces.expr,1,function(x,all,selected){
	yes <- x[all%in%selected]
	no <- x[!all%in%selected]
	test <- t.test(yes,no)$p.value
	},all=colnames(puces.expr),selected=row.names(df.LMS.extract)[df.LMS.extract$Cluster=="hLMS"])
ttest.adj <- p.adjust(ttest.test)
ttest.sig <- -log10(ttest.adj)

t.stable <- apply(puces.expr,1,function(x,df.LMS){
	stable <- x[df.LMS$Cluster=="hLMS"]
	other <- x[df.LMS$Cluster=="oLMS"]
	score <- (mean(stable)-mean(other)) / sqrt((sd(stable)^2)/length(stable) + (sd(other)^2)/length(other))
	},df.LMS=df.LMS.extract)

var.stable <- apply(puces.expr[,row.names(df.LMS.extract)[df.LMS.extract$Cluster=="hLMS"]],1,var)
var.other <- apply(puces.expr[,row.names(df.LMS.extract)[df.LMS.extract$Cluster=="oLMS"]],1,var)

median.stable <- apply(puces.expr[,row.names(df.LMS.extract)[df.LMS.extract$Cluster=="hLMS"]],1,median)
median.other <- apply(puces.expr[,row.names(df.LMS.extract)[df.LMS.extract$Cluster=="oLMS"]],1,median)


df_tests_diff <- data.frame(ttest.adj=ttest.sig,FDR=-log10(p.adjust(ttest.test,"fdr")),logRatio=t.stable,row.names=names(t.stable),median.hLMS=median.stable,median.oLMS=median.other,var.hLMS=var.stable,var.oLMS=var.other)
df_tests_diff <- df_tests_diff[order(df_tests_diff$logRatio,decreasing=T),]

write.table(data.frame(SYMBOL=row.names(df_tests_diff),t=df_tests_diff$logRatio),file.path(res.dir,"Affy_t_scores_OCEANCODE.rnk"),sep="\t",row.names=F,col.names=F,quote=F)
writeLines(row.names(df_tests_diff)[df_tests_diff$ttest.adj>2 & df_tests_diff$logRatio>0],file.path(res.dir,"hLMS_genes_up.txt"))
writeLines(row.names(df_tests_diff)[df_tests_diff$ttest.adj>2 & df_tests_diff$logRatio<0],file.path(res.dir,"hLMS_genes_down.txt"))

median.rows <- apply(puces.expr,1,median)
sd.rows <- apply(puces.expr,1,sd)
dataset.center <- (puces.expr - median.rows)/sd.rows

mat <- puces.expr[row.names(puces.expr)%in%row.names(df_tests_diff[df_tests_diff$ttest.adj>2,]),]
res.pca <- PCA(t(mat),scale.unit=FALSE,ncp=3, graph = FALSE)
res.hcpc <- HCPC(res.pca, nb.clust=2, conso=0,graph = FALSE)

mat.temp <- dataset.center[,match(res.hcpc$call$t$tree$labels,colnames(dataset.center))]

#mat.temp <- dataset.center
clust <- pheatmap(mat.temp[row.names(mat.temp)%in%row.names(df_tests_diff[df_tests_diff$ttest.adj>2,]),],annotation_col=df.LMS.extract,cluster_rows=T, cluster_cols=res.hcpc$call$t$tree,
	annotation_colors=generate.annot.col(df.LMS.extract),show_rownames = F,show_colnames = F,#clustering_distance_cols="correlation",clustering_method="ward",
	color=colorRampPalette(c(colours()[128],colours()[128],"white","red","red"))(100),breaks=seq(-6,6,0.12),silent=T)
cluster=factor(cutree(clust$tree_col, k = 2))
df.LMS.extract$Cluster2 <- cluster[match(row.names(df.LMS.extract),names(cluster))]
df.LMS.extract$Cluster2 <- ifelse(df.LMS.extract$Cluster2==2,"hLMS","oLMS")
df.LMS <- df.LMS.extract
df.LMS$Cluster <- df.LMS.extract$Cluster2
df.LMS$Cluster2 <- NULL
df.LMS$Type <- NULL
colnames(df.LMS) <- sub("Localisation","Location",colnames(df.LMS))

pdf(file.path(res.dir,"Figure1B_DE_Affymetrix.pdf"), paper="special",width=9,height=7)
pheatmap(mat.temp[row.names(mat.temp)%in%row.names(df_tests_diff[df_tests_diff$ttest.adj>2,]),],annotation_col=df.LMS,cluster_rows=T, cluster_cols=res.hcpc$call$t$tree,
	annotation_colors=generate.annot.col(df.LMS),show_rownames = F,show_colnames = F,#clustering_distance_cols="correlation",clustering_method="ward",
	color=colorRampPalette(c(colours()[128],colours()[128],"white","red","red"))(100),breaks=seq(-6,6,0.12),cutree_col=2,)
#export.plot(file.path(res.dir,"Figure1B_DE_Affymetrix"),width=9,heigh=7)
dev.off()
message("Figure 1B OK")

df.LMS <- clinical
m <- match(row.names(df.LMS),row.names(df.LMS.extract))
df.LMS <- df.LMS[!is.na(m),]
df.LMS$Cluster <- df.LMS.extract$Cluster2[na.omit(m)]
save(df.LMS,df_tests_diff,puces.expr,file=file.path(res.dir,"DE_LMS_affy_objects.RData"))

m <- match(row.names(df.cols),row.names(df.LMS))
df.cols$Cluster <- df.LMS$Cluster[m]

	pdf(file.path(res.dir,"Figure1A_all_samples_clustering.pdf"), paper="special",width=11,height=15)
pheatmap(harmonized.scaled,annotation_col=df.cols[,-c(7)],annotation_row=df.row,
annotation_colors=c(generate.annot.col(df.cols[,-c(7)]),generate.annot.col(df.row)),cluster_rows=F,cluster_cols=res.hcpc.harmonized$call$t$tree,
	show_rownames = F,show_colnames = F,color=colorRampPalette(c(colours()[128],colours()[128],"white","red","red"))(100),breaks=seq(-7,7,0.14))
#export.plot(file.path(res.dir,"Figure1A_all_samples_clustering"),width=11,heigh=15)
dev.off()
message("Figure 1A OK")







