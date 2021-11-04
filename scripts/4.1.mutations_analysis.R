source("4.functions.R")

message("Start distribution mutation analysis")
cat <- read.table(file.path(data.dir,"ext","TP53_RB1_status.tab"),sep="\t",head=T)
cat$Code_fragment <- sub("T$","",as.vector(cat$Code_fragment))
cat$Code_fragment <- sub("T","_",as.vector(cat$Code_fragment))
cat$Cluster <- merge.clinical$Cluster[match(cat$Code_fragment,row.names(merge.clinical))]

cat.grouped <- cat
cat.grouped$TP53 <- NA
cat.grouped[cat$TP53%in%c("MS/FS","MS/L","FS/L"),]$TP53 <- "2T"
cat.grouped[cat$TP53%in%c("2L","MS"),]$TP53 <- "1T"

cat.grouped$RB1 <- NA
cat.grouped[cat$RB1%in%c("MUT/L"),]$RB1 <- "2T"
cat.grouped[cat$RB1%in%c("2L","MUT"),]$RB1 <- "1T"

message("Fisher exact test on TP53 and RB1 grouped effect by type")
test.cat.grouped <- rbindlist(lapply(colnames(cat)[2:3],clinical.fisher,df.LMS=cat.grouped))
print(test.cat.grouped)

message("Fisher exact test on TP53 and RB1")
test.cat <- rbindlist(lapply(colnames(cat)[2:3],clinical.fisher,df.LMS=cat))
print(test.cat)

cat.plot <- rbindlist(lapply(colnames(cat)[2:3],plot.prop.mutation,df=cat))

cat.plot$gene <- factor(as.vector(cat.plot$gene),levels=c("TP53","RB1"))
g <- ggplot(cat.plot[cat.plot$gene%in%c("TP53","RB1")],aes(x=alteration,y=Cluster)) + 
		geom_point(aes(size=sample.freq,color=Cluster,fill=Cluster),shape=21,alpha=0.5) +
		scale_fill_manual(values=c(hLMS=colours()[613],oLMS=colours()[131])) +
		scale_color_manual(values=c(hLMS=colours()[613],oLMS=colours()[131])) +
		theme_pubr() + scale_size(range=c(1,10),breaks=seq(0,1,0.1)) + ylab("") + facet_wrap(~gene,ncol=2,scales="free_x")

pdf(file.path(res.dir,"Figure4C_cat_alteration_categories.pdf"),width=5,heigh=3)
print(g)
#export.plot(file.path(res.dir,paste0("Figure4C_cat_alteration_categories")),width=5,height=3)
dev.off()
message("Figure 4C OK")



valid <- read.table(file.path(data.dir,"ext","Table_mut_6_genes.tab"),sep="\t",head=T)
valid.loc <- valid[,c("DMD.All","DMD.Dp427","PTEN.Loc","ATRX.Loc","Cluster")]

message("Fisher exact test on ATRX, PTEN, DMD cellular localisation")
test.loc <- rbindlist(lapply(colnames(valid.loc)[1:4],clinical.fisher,df.LMS=valid.loc))
print(test.loc)


valid.plot <- rbindlist(lapply(colnames(valid.loc)[-5],plot.prop.mutation,df=valid.loc))
valid.plot$gene <- factor(as.vector(valid.plot$gene),levels=colnames(valid.loc)[-5])
valid.plot <- valid.plot[!is.na(valid.plot$alteration),]
g <- ggplot(valid.plot,aes(x=alteration,y=Cluster)) + 
		geom_point(aes(size=sample.freq,color=Cluster,fill=Cluster),shape=21,alpha=0.5) +
		scale_fill_manual(values=c(hLMS=colours()[613],oLMS=colours()[131])) +
		scale_color_manual(values=c(hLMS=colours()[613],oLMS=colours()[131])) + xlab("signal") +
		theme_pubr() + scale_size(range=c(1,10),breaks=seq(0,1,0.1)) + ylab("") + facet_wrap(~gene,ncol=2,scales="free_x")#,scales="free")
pdf(file.path(res.dir,"Figure4E_loc_alteration_categories.pdf"),width=5,heigh=6)
print(g)
#export.plot(file.path(res.dir,paste0("Figure4E_loc_alteration_categories")),width=5,height=6)
dev.off()
message("Figure 4E OK")

#####################################################################
message("Analysis of DMD isoform expression in ICGC cohort")

isoform_DMD <- read.table(file.path(data.dir,"ext","DMD_isoform_expression_ICGC.tab"),sep="\t",head=T,row=1)

isoform_DMD$Cluster <- merge.clinical$Cluster[match(row.names(isoform_DMD),row.names(merge.clinical))]


isoform_DMD <- melt(isoform_DMD)
isoform_DMD <- isoform_DMD[isoform_DMD$variable%in%c("Dp427","Dp45","Dp71"),]
isoform_DMD$variable <- factor(as.vector(isoform_DMD$variable),levels=c("Dp427","Dp45","Dp71"))

g <- ggplot(isoform_DMD,aes(x=Cluster,y=value)) + geom_violin(aes(color=Cluster)) +
	theme_pubr() + scale_color_manual(values=c(hLMS=colours()[613],oLMS=colours()[131])) + geom_jitter(width=0.1,aes(color=Cluster)) +
	facet_wrap(~variable,scales="free",ncol=3) + stat_compare_means(comparison=list(c("hLMS","oLMS")),method="t.test") 
pdf(file.path(res.dir,"Figure4D_DMD_isoforms.pdf"),width=5.5,heigh=2.5)
print(g)
#export.plot(file.path(res.dir,"Figure4D_DMD_isoforms"),width=5.5,heigh=2.5)
dev.off()
message("Figure 4D OK")


#####################################################################
message("Start analysis mutational patterns")

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
load(file.path(data.dir,"RData","vcfs_ICGC_all_somatic.RData")) # vcfs

message("Compute 96 trinucleotide patterns")
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)

mut_mat <- mut_mat + 0.0001

message("Retrieve COSMIC signature from https://cancer.sanger.ac.uk/cancergenome/assets/")
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/","signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures = cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
signatures30 = as.matrix(cancer_signatures[,4:33])

message("Compute cosine similarity")
cos_sim_samples_signatures = cos_sim_matrix(mut_mat, signatures30)
cos_sim_ss <- cos_sim_samples_signatures[,apply(cos_sim_samples_signatures,2,max)>0.75]

hclust_cosmic = cluster_signatures(signatures30[,colnames(cos_sim_ss)], method = "average")
cosmic_order = colnames(signatures30)[hclust_cosmic$order]

pdf(file.path(res.dir,"FigureS4E_cosine_cosmic_signature.pdf"),width=6,heigh=7)
clust <- pheatmap(cos_sim_ss,
	annotation_row=data.frame(Cluster=annots.ICGC$Cluster,row.names=row.names(annots.ICGC)),
	clustering_method="complete",cluster_cols=hclust_cosmic,
	breaks=seq(0.74,1,0.0026),color=colorRampPalette(c("white","red","firebrick","brown","black"))(100),
	annotation_colors=generate.annot.col(data.frame(Cluster=annots.ICGC$Cluster,row.names=row.names(annots.ICGC))))
#export.plot(file.path(res.dir,"FigureS4E_cosine_cosmic_signature"),width=6,heigh=7)
dev.off()
message("Figure S4E heatmap OK")



message("Compute contribution of COSMIC signatures")
fit_res <- fit_to_signatures(as.matrix(mut_mat), as.matrix(signatures30))
cosmic_contrib <- t(fit_res$contribution)/rowSums(t(fit_res$contribution))
cosmic_contrib <- cosmic_contrib[,colnames(cosmic_contrib)%in%colnames(cos_sim_ss)]

cosmic_contrib.m <- as.data.frame(cosmic_contrib)
cosmic_contrib.m$Patients <- row.names(cosmic_contrib.m)
cosmic_contrib.m <- melt(cosmic_contrib.m)
colnames(cosmic_contrib.m)[3] <- "contribution"

cos_sim_ss.m <- as.data.frame(cos_sim_ss)
cos_sim_ss.m$Patients <- row.names(cos_sim_ss.m)
cos_sim_ss.m <- melt(cos_sim_ss.m)
colnames(cos_sim_ss.m)[3] <- "cosine"

test <- merge(cos_sim_ss.m,cosmic_contrib.m,by=c("Patients","variable"))
test$Cluster <- annots.ICGC$Cluster[match(test$Patients,row.names(annots.ICGC))]
test <- test[test$cosine>0.75,]
hclust_contrib = hclust(dist(cos_sim_ss),method="complete",)

test$Patients <- factor(as.vector(test$Patients),levels=row.names(cosmic_contrib)[clust$tree_row$order])
g <- ggplot(test,aes(x=Patients,y=contribution,group=Patients)) + 
	geom_col(aes(fill=variable),width=0.8) + scale_fill_manual(values=brewer.pal(8, "Set3")) +
	coord_flip() + theme(panel.background = element_rect(fill = "white", colour = "black")) + ylim(c(0,1))
pdf(file.path(res.dir,"FigureS4E_contribution_cosmic_signature_barplot.pdf"),width=4,heigh=7)
print(g)
#export.plot(file.path(res.dir,"FigureS4E_contribution_cosmic_signature_barplot"),width=4,heigh=7)
dev.off()
message("Figure S4E barplot OK")

###########################################################################

message("Compute tumor mutational burden")
nbMB <- sum(seqlengths(get(ref_genome))[1:24])
TMB <- lapply(vcfs,function(x,nbMB){ (length(x)/nbMB)*10^6},nbMB)
m <- match(row.names(annots.ICGC),names(TMB))
annots.ICGC$TBM <- unlist(TMB)[m]


message("Load breakpoint prediction")
load(file.path(data.dir,"RData","fusion_breakpoints.RData"))
bp.coords <- bp.coords[names(bp.coords)%in%row.names(merge.clinical)]
counter <- 1
bp.counts <- rbindlist(lapply(bp.coords,function(x,samp){
	res <- data.table(Patient=samp[counter],nb.bp=length(x$pos1))
	counter <<- counter + 1
	res
	},names(bp.coords)))


m <- match(row.names(annots.ICGC),bp.counts$Patient)
annots.ICGC$BP <- bp.counts$nb.bp[m]

mut <- melt(annots.ICGC[row.names(annots.ICGC)!="LMS23",c("TBM","BP","Cluster")])

g <- ggplot(mut,aes(x=Cluster,y=value)) + geom_boxplot(aes(color=Cluster)) + facet_wrap(~variable,scales="free_y",ncol=1) +
	scale_color_manual(values=c(hLMS=colours()[613],oLMS=colours()[131])) +
	theme_bw() + stat_compare_means(comparisons=list(c("hLMS","oLMS")))
pdf(file.path(res.dir,"FigureS4D_TBM_BP.pdf"),width=3,heigh=6)
print(g)
#export.plot(file.path(res.dir,"FigureS4D_TBM_BP"),width=3,heigh=6)
dev.off()
message("Figure S4D boxplots OK")








