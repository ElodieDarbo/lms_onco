message("Start miRComb miRNA-mRNA analysis")

message("Load databases")
load(file.path(data.dir,"RData","all_db_miRNA_target.RData"))
miRecords=miRNA_targets$miRecords
miRTarBase=miRNA_targets$miRTarBase
miRecords$mir17 <- gsub("\\[|\\]","",as.vector(miRecords$mir17))
miRecords$target_name <- gsub("-","",miRecords$target_name)
##############################################################
#### all diff
diff.miRNA <- miRNA_ICGC_TCGA[Significance=="TCGA+ICGC"]$ID
hairpins <- tolower(sub("-[35]p","",diff.miRNA))

hairpins.TCGA <- tolower(sub("-[35]p","",row.names(miRNA_TCGA)))
hairpins.ICGC <- tolower(sub("-[35]p","",row.names(miRNA_ICGC_norm)))

message("Compute average expression per pre-miRNA")
message("ICGC")

miRNA_ICGC_norm <- miRNA_ICGC_norm[,match(row.names(annots.ICGC),colnames(miRNA_ICGC_norm))]
miRNA_ICGC <- miRNA_ICGC_norm[row.names(miRNA_ICGC_norm)%in%diff.miRNA,]
miRNA_ICGC.hairpins <- as.data.frame(miRNA_ICGC)
hairpins.ICGC <- tolower(sub("-[35]p","",row.names(miRNA_ICGC)))
miRNA_ICGC.hairpins$hairpins <- hairpins.ICGC
miRNA_ICGC.hairpins <- as.data.frame(data.table(miRNA_ICGC.hairpins)[,lapply(.SD,mean),by=hairpins])
row.names(miRNA_ICGC.hairpins) <- miRNA_ICGC.hairpins$hairpins
miRNA_ICGC.hairpins <- miRNA_ICGC.hairpins[,-1]

all.equal(row.names(annots.ICGC),colnames(ICGC.expr))
all.equal(row.names(annots.ICGC),colnames(miRNA_ICGC.hairpins))
all.equal(colnames(miRNA_ICGC.hairpins),colnames(ICGC.expr))

message("TCGA")
miRNA_TCGA <- miRNA_TCGA[,colnames(miRNA_TCGA)%in%row.names(annots.TCGA) & colnames(miRNA_TCGA)%in%colnames(TCGA.expr)]
annots.TCGA <- annots.TCGA[row.names(annots.TCGA)%in%colnames(miRNA_TCGA),]
m <- match(row.names(annots.TCGA),colnames(miRNA_TCGA))
miRNA_TCGA <- miRNA_TCGA[,m]
m <- match(row.names(annots.TCGA),colnames(TCGA.expr))
TCGA.expr <- TCGA.expr[,m]

miRNA_TCGA <- miRNA_TCGA[row.names(miRNA_TCGA)%in%diff.miRNA,]
hairpins.TCGA <- tolower(sub("-[35]p","",row.names(miRNA_TCGA)))
miRNA_TCGA.hairpins <- as.data.frame(miRNA_TCGA)
miRNA_TCGA.hairpins$hairpins <- hairpins.TCGA
miRNA_TCGA.hairpins <- as.data.frame(data.table(miRNA_TCGA.hairpins)[,lapply(.SD,mean),by=hairpins])
row.names(miRNA_TCGA.hairpins) <- miRNA_TCGA.hairpins$hairpins
miRNA_TCGA.hairpins <- miRNA_TCGA.hairpins[,-1]

miRecords$mir17 <- tolower(sub("-[35]p","",miRecords$mir17))
miRecords <- unique(miRecords)
miRTarBase <- unique(miRTarBase)
row.names(miRTarBase) <- miRTarBase$names


message("Run miRComb")
message("ICGC")
ICGC.expr[is.na(ICGC.expr)] <- 0
miRAll.ICGC <- run_mirComb(annots.ICGC,ICGC.expr,miRNA_ICGC.hairpins,diff.miRNA,DIO3=as.vector(miRNA_ICGC_TCGA[DIO3_DLK1=="yes"]$ID))
message("TCGA")
miRAll.TCGA <- run_mirComb(annots.TCGA,TCGA.expr,miRNA_TCGA.hairpins,diff.miRNA,DIO3=as.vector(miRNA_ICGC_TCGA[DIO3_DLK1=="yes"]$ID))

net.ICGC <- data.table(miRAll.ICGC@net)[adj.pval<0.01 ]
net.TCGA <- data.table(miRAll.TCGA@net)[adj.pval<0.01 ]

net.ICGC[dat.sum>=0,ID:=paste(miRNA,mRNA,sep=":")]
net.TCGA[dat.sum>=0,ID:=paste(miRNA,mRNA,sep=":")]

net.merged <- merge(net.ICGC[,list(ID,miRNA,mRNA,dat.miRecords,dat.miRTarBase,dat.sum,logFC.mRNA.ICGC=logratio.mRNA,logFC.miRNA.ICGC=logratio.miRNA,cor.ICGC=cor)],
                      net.TCGA[,list(ID,miRNA,mRNA,dat.miRecords,dat.miRTarBase,dat.sum,logFC.mRNA.TCGA=logratio.mRNA,logFC.miRNA.TCGA=logratio.miRNA,cor.TCGA=cor)],
                      by=c("ID","miRNA","mRNA","dat.miRecords","dat.miRTarBase","dat.sum"))
net.merged$DIO3_DLK1 <- ifelse(net.merged$miRNA%in%DIO3,"yes","no")
write.table(net.merged,file.path(res.dir,"miRComb_table.tab"),sep="\t",row.names=F,quote=F)
message("Table S2C OK")

stats.genes <- net.merged[dat.sum>0,list(nb.miRNA=length(miRNA),nb.DIO3=length(miRNA[DIO3_DLK1=="yes"]),sign=ifelse(unique(logFC.mRNA.ICGC)>0,"up","dn")),by=c("mRNA")]

writeLines(as.vector(stats.genes$mRNA),file.path(res.dir,"test.mRNA.miRNA.txt"))

message("import results from GSAn and visualize")
GSAN <- read.table(file.path(data.dir,"precomp","GSAn_result_mirDiff.csv"),sep="\t",head=T,quote="\"")
GSAN <- GSAN[GSAN$Term_depth>5,]
GSAN.members <- rbindlist(apply(GSAN,1,extract.genes))
GSAN.members$group <- ifelse(GSAN.members$SYMBOL%in%stats.genes[sign=="up"]$mRNA,"up","down")
GSAN.members$DIO3 <- ifelse(GSAN.members$SYMBOL%in%stats.genes[nb.DIO3>0]$mRNA,"yes","no")
GSAN.members$no.DIO3 <- ifelse(GSAN.members$SYMBOL%in%stats.genes[nb.miRNA>nb.DIO3]$mRNA,"yes","no")

#GSAN.members <- GSAN.members[,list(nb.up.DIO3=length(SYMBOL[no.DIO3=="yes" & DIO3=="yes"])/nrow(stats.genes[sign=="up"]),nb.up=length(SYMBOL[group=="up" & DIO3!="yes"])/nrow(stats.genes[sign=="up"]),nb.DIO3=length(SYMBOL[DIO3=="yes" & no.DIO3=="no"])/nrow(stats.genes[sign=="up"]),nb.down=-length(SYMBOL[group=="down"])/nrow(stats.genes[sign=="dn"])),by=c("term","ont")]
GSAN.members <- GSAN.members[,list(nb.up.DIO3=length(SYMBOL[no.DIO3=="yes" & DIO3=="yes"]),nb.up=length(SYMBOL[group=="up" & DIO3!="yes"]),nb.DIO3=length(SYMBOL[DIO3=="yes" & no.DIO3=="no"]),nb.down=-length(SYMBOL[group=="down"])),by=c("term","ont")]

GSAN.members <- melt(GSAN.members)
GSAN.members$variable <- factor(as.vector(GSAN.members$variable),levels=c("nb.up","nb.up.DIO3","nb.DIO3","nb.down"))

GSAN.members$term <- sub("RNA polymerase II","RNAPolII",GSAN.members$term)
GSAN.members <- data.table(GSAN.members)
GSAN.members$term <- factor(as.vector(GSAN.members$term),levels=unique(GSAN.members[variable!="nb.down"][,list(value=sum(value),ont),by=term][order(value)][order(ont)]$term))


GSAN.members$grouped.term <- "aother"
GSAN.members[term%in%c("plasma membrane bounded cell projection assembly","extracellular matrix disassembly","branching morphogenesis of an epithelial tube")]$grouped.term <- "zcell.migration"
GSAN.members[term%in%c("regulation of heart contraction","cation channel activity","calcium ion transport","inorganic cation transmembrane transport")]$grouped.term <- "ycell.contraction"
GSAN.members[term%in%c("cell cycle arrest")]$grouped.term <- "ycell.cycle"
GSAN.members[term%in%c("RNAPolII cis-regulatory region sequence-specific DNA binding","transcription by RNAPolII")]$grouped.term <- "xtranscriptiolan.regulation"
GSAN.members[term%in%c("extracellular exosome","cytoplasmic vesicle","secretory granule membrane","lysosomal membrane","neutrophil degranulation")]$grouped.term <- "tcellular.transport"
GSAN.members[term%in%c("ATP binding","protein serine/threonine kinase activity","cellular protein modification process")]$grouped.term <- "oprotein.modification"

GSAN.members$grouped.color <- "black"
GSAN.members[term%in%c("ATP binding","protein serine/threonine kinase activity","cellular protein modification process")]$grouped.color <-  colours()[72]
GSAN.members[term%in%c("extracellular exosome","cytoplasmic vesicle","secretory granule membrane","lysosomal membrane","neutrophil degranulation")]$grouped.color <-  "black"
GSAN.members[term%in%c("RNAPolII cis-regulatory region sequence-specific DNA binding","transcription by RNAPolII")]$grouped.color <- colours()[72]
GSAN.members[term%in%c("regulation of heart contraction","cation channel activity","calcium ion transport","inorganic cation transmembrane transport")]$grouped.color <- "black"
GSAN.members[term%in%c("cell cycle arrest")]$grouped.color <- colours()[72]
GSAN.members[term%in%c("plasma membrane bounded cell projection assembly","extracellular matrix disassembly","branching morphogenesis of an epithelial tube")]$grouped.color <- "black"


GSAN.members$term <- factor(as.vector(GSAN.members$term),levels=unique(GSAN.members[order(grouped.term),]$term))
GSAN.members <- GSAN.members[order(grouped.term)]

grouped.color <- GSAN.members[,list(coloring=unique(grouped.color)),by=term]$coloring

g <- ggplot(GSAN.members[term!="nucleus"],aes(y=value,x=term)) + geom_col(aes(fill=variable)) + 
  coord_flip() + ylim(-30,40)  + scale_fill_manual(values=c(nb.up=colours()[34],nb.up.DIO3=colours()[35],nb.DIO3=colours()[36],nb.down=colours()[121])) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),axis.text.y = element_text(colour = grouped.color) ) 


pdf(file.path(res.dir,"Figure2D_GSAN_miRNA_targets.pdf"),width=8,heigh=4)
print(g)
#export.plot(file.path(res.dir,"Figure2D_GSAN_miRNA_targets"),width=8,heigh=4)
dev.off()
message("Figure 2D OK")

