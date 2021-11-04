copy.test <- function(x){
  x <- as.character(x)
  y <- x
  y <- sub("A",4,y)
  y <- sub("G",3,y)
  y <- sub("N",2,y)
  y <- sub("P",1,y)
  y <- sub("D",0,y)
  y <- sub("L",1,y)
  multi <- grep("-",y)
  y[multi] <- lapply(multi,function(m,x){
    tmp <- min(as.numeric(unlist(strsplit(x[multi],"-"))))

    },y)
  y <- as.numeric(y)
  return(y)
}

import.cyto <- function(data.dir,Txdb,org.Hs.eg.db){
  cytoband.hg19 <- read.table(file.path(data.dir,"cytoBand_hg19.txt"),sep="\t")
  all.hg19 <- genes(Txdb)
  all.hg19 <- all.hg19[order(start(all.hg19),decreasing=F)]
  band <- sapply(strsplit(as.vector(cytoband.hg19$V4),"[.]"),"[[",1)
  band <- as.numeric(sub("[pq]","",band))
  arm <- sapply(strsplit(as.vector(cytoband.hg19$V4),""),"[[",1)
  cytoband.hg19 <- data.table(chromosome=cytoband.hg19$V1,start=cytoband.hg19$V2,end=cytoband.hg19$V3,arm=arm,band=band)
  cytoband.hg19 <- cytoband.hg19[,list(start=min(start),end=max(end)),by=c("chromosome","arm","band")]
  cytoband.hg19 <- GRanges(cytoband.hg19$chromosome,IRanges(cytoband.hg19$start,cytoband.hg19$end),arm=cytoband.hg19$arm,band=cytoband.hg19$band)
  geneIDs <- bitr(all.hg19$gene_id,fromType="ENTREZID",toType="SYMBOL",org.Hs.eg.db)
  m <- match(all.hg19$gene_id,geneIDs$ENTREZID)
  all.hg19$SYMBOL <- geneIDs$SYMBOL[m]
  ovl <- findOverlaps(cytoband.hg19,all.hg19)
  cytoband.hg19 <- cytoband.hg19[queryHits(ovl)]
  cytoband.hg19$SYMBOL <- all.hg19$SYMBOL[subjectHits(ovl)]
  cytoband.hg19$gene.start <- start(all.hg19)[subjectHits(ovl)]
  cytoband.hg19$gene.end <- end(all.hg19)[subjectHits(ovl)]
  cytoband.hg19 <- cytoband.hg19[!is.na(cytoband.hg19$SYMBOL)]
  cytoband.hg19 <- cytoband.hg19[!duplicated(cytoband.hg19$SYMBOL)]
  cytoband.hg19$ont <- paste(tolower(as.vector(seqnames(cytoband.hg19))),cytoband.hg19$arm,cytoband.hg19$band,sep="")
  cytoband.hg19$chrNum <- sub("chr","",as.vector(seqnames(cytoband.hg19)))
  cytoband.hg19$chrNum[cytoband.hg19$chrNum=="X"] <- 23
  cytoband.hg19$chrNum[cytoband.hg19$chrNum=="Y"] <- 24
  cytoband.hg19$chrNum <- as.numeric(cytoband.hg19$chrNum)
  return(cytoband.hg19)
}

enrich.CN <- function(x,classif,group=NULL){
  classif <- as.vector(classif) 
  if (is.null(group)){
    classif <- classif[!is.na(x)]
    x <- x[!is.na(x)]
    LMS <- which(classif=="selected")
    other <- which(classif!="selected")
    loss <- which(x==1)
  gain <- which(x==3)
  stable <- which(x==2)
  amp <- which(x==4)
  del <- which(x==0)
  mat <- matrix(c(length(which(LMS%in%loss)),length(which(!LMS%in%loss)),length(which(other%in%loss)),length(which(!other%in%loss))),ncol=2,byrow = T)
  fisher.loss <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(length(which(LMS%in%gain)),length(which(!LMS%in%gain)),length(which(other%in%gain)),length(which(!other%in%gain))),ncol=2,byrow = T)
  fisher.gain <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(length(which(LMS%in%stable)),length(which(!LMS%in%stable)),length(which(other%in%stable)),length(which(!other%in%stable))),ncol=2,byrow = T)
  fisher.stable <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(length(which(LMS%in%del)),length(which(!LMS%in%del)),length(which(other%in%del)),length(which(!other%in%del))),ncol=2,byrow = T)
  fisher.del <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(length(which(LMS%in%amp)),length(which(!LMS%in%amp)),length(which(other%in%amp)),length(which(!other%in%amp))),ncol=2,byrow = T)
  fisher.amp <- fisher.test(mat,alternative = "greater")
  
  }
  else{
    l <- length(classif)
    del <- x[1:l]
    loss <- x[(l+1):(2*l)]
    stable <- x[(l*2+1):(3*l)]
    gain <- x[(l*3+1):(4*l)]
    amp <- x[(l*4+1):(5*l)]
    idx <- classif=="selected"
    mat <- matrix(c(sum(loss[idx]),sum(c(stable[idx],gain[idx],del[idx],amp[idx])),
      sum(loss[!idx]),sum(c(stable[!idx],gain[!idx],del[!idx],amp[!idx]))),ncol=2,byrow = T)
  fisher.loss <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(sum(gain[idx]),sum(c(stable[idx],loss[idx],del[idx],amp[idx])),
      sum(gain[!idx]),sum(c(stable[!idx],loss[!idx],del[!idx],amp[!idx]))),ncol=2,byrow = T)
  fisher.gain <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(sum(stable[idx]),sum(c(loss[idx],gain[idx],del[idx],amp[idx])),
      sum(stable[!idx]),sum(c(loss[!idx],gain[!idx],del[!idx],amp[!idx]))),ncol=2,byrow = T)
  fisher.stable <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(sum(del[idx]),sum(c(loss[idx],gain[idx],stable[idx],amp[idx])),
      sum(del[!idx]),sum(c(loss[!idx],gain[!idx],stable[!idx],amp[!idx]))),ncol=2,byrow = T)
  fisher.del <- fisher.test(mat,alternative = "greater")
  mat <- matrix(c(sum(amp[idx]),sum(c(loss[idx],gain[idx],stable[idx],del[idx])),
      sum(amp[!idx]),sum(c(loss[!idx],gain[!idx],stable[!idx],del[!idx]))),ncol=2,byrow = T)
  fisher.amp <- fisher.test(mat,alternative = "greater")

  }
  res <- c(fisher.del$p.value,fisher.loss$p.value,fisher.stable$p.value,fisher.gain$p.value,fisher.amp$p.value)
  return(res)
}


import.copyNumber.by.genes <- function(CN.puces,df.CN,cytoband.hg19,logRatio,thr.z,fisher=TRUE){
  message("start import Copy Number data")
  genes.puces <- data.frame(ont=cytoband.hg19$ont,gene=cytoband.hg19$SYMBOL,chr=as.vector(seqnames(cytoband.hg19)),arm=cytoband.hg19$arm,band=cytoband.hg19$band,chrNum=cytoband.hg19$chrNum)
  CN.puces <- CN.puces[row.names(CN.puces)%in%genes.puces$gene,]
  genes.puces <- genes.puces[genes.puces$gene%in%row.names(CN.puces),]
  CN.puces <- CN.puces[match(as.vector(genes.puces$gene),row.names(CN.puces)),]
  selected.samples <- row.names(df.CN)[df.CN$Cluster=="hLMS"]
  other.samples <- row.names(df.CN)[df.CN$Cluster!="hLMS"]
  CN.stable <- CN.puces[,selected.samples]
  CN.other <- CN.puces[,other.samples]
  genes.puces$lost.stable <- apply(CN.stable,1,function(x){
  sum(!is.na(x) & x < 2)
  })
  genes.puces$gain.stable <- apply(CN.stable,1,function(x){
  sum(!is.na(x) & x > 2)
  })
  genes.puces$lost.other <- apply(CN.other,1,function(x){
  sum(!is.na(x) & x < 2)
  })
  genes.puces$gain.other <- apply(CN.other,1,function(x){
  sum(!is.na(x) & x > 2)
  })
  genes.puces$del.stable <- apply(CN.stable,1,function(x){
  sum(!is.na(x) & x == 0)
  })
  genes.puces$amp.stable <- apply(CN.stable,1,function(x){
  sum(!is.na(x) & x == 4)
  })
  genes.puces$del.other <- apply(CN.other,1,function(x){
  sum(!is.na(x) & x ==0 )
  })
  genes.puces$amp.other <- apply(CN.other,1,function(x){
  sum(!is.na(x) & x == 4)
  })
  genes.puces$noNa.other <- apply(CN.other,1,function(x){
  sum(!is.na(x))
  })
  genes.puces$noNa.stable <- apply(CN.stable,1,function(x){
  sum(!is.na(x))
  })
  genes.puces.dt <- data.table(genes.puces)
  CN.puces <- CN.puces[!duplicated(genes.puces.dt$gene),]
  genes.puces.dt <- genes.puces.dt[!duplicated(genes.puces.dt$gene)]
  m <- match(genes.puces.dt$gene,names(logRatio))
  genes.puces.dt$z.expr <- logRatio[m]
  if (fisher){
    annots.genes <- as.data.frame(genes.puces.dt[,list(arm,chrNum=factor(chrNum),z.expr)])
    row.names(annots.genes) <- genes.puces.dt$gene
    chrom.gaps <- cumsum(Rle(genes.puces.dt$chrNum)@lengths)
    message("start Fisher test genes")
    testfisher <- t(apply(CN.puces[,row.names(df.CN)],1,enrich.CN,classif=ifelse(df.CN$Cluster=="hLMS","selected","other")))
    testfisher <- -log(testfisher,base=10)
    testfisher <- as.data.frame(testfisher)
    colnames(testfisher) <- c("deletion","loss","stable","gain","amplification")
    testfisher <- as.data.frame(testfisher)
    testfisher$z.expr <- genes.puces.dt$z.expr
    genes.puces.dt$loss.gene <- testfisher$loss
    genes.puces.dt$gain.gene <- testfisher$gain
    genes.puces.dt$del.gene <- testfisher$deletion
    genes.puces.dt$amp.gene <- testfisher$amplification
    thr <- -log10(thr.z)
    annots.genes$fisher_genes <- 0
    annots.genes$fisher_genes[testfisher$loss> thr] <- -1
    annots.genes$fisher_genes[testfisher$gain> thr] <- 1
    annots.genes$fisher_genes[testfisher$del> thr & apply(cbind(testfisher$loss,testfisher$gain,testfisher$amp),1,max) < testfisher$del] <- -2
    annots.genes$fisher_genes[testfisher$amp> thr & apply(cbind(testfisher$loss,testfisher$gain,testfisher$del),1,max) < testfisher$amp] <- 2 
    annots.genes$fisher_genes <- factor(annots.genes$fisher_genes,levels=c(-2,-1,0,1,2))
    annots.genes$chrNum <- factor(annots.genes$chrNum,levels=1:24)
    chrom.gaps <- cumsum(Rle(sort(as.numeric(as.vector(annots.genes$chrNum))))@lengths)
    res <- list(CN=CN.puces,annots.genes=annots.genes,chrom.gaps=chrom.gaps,genes.dt=genes.puces.dt,annots.samples=df.CN)
  } else {
    res <- list(genes.dt=genes.puces.dt)
  }
  return(res)
}

percentage.events <- function(genes.dt.all,cytoband.hg19){
  genes.dt.all$percent.gain.h <- (genes.dt.all$gain.stable-genes.dt.all$amp.stable)/genes.dt.all$noNa.stable
  genes.dt.all$percent.amp.h <- genes.dt.all$amp.stable/genes.dt.all$noNa.stable
  genes.dt.all$percent.lost.h <- (genes.dt.all$lost.stable-genes.dt.all$del.stable)/genes.dt.all$noNa.stable
  genes.dt.all$percent.del.h <- genes.dt.all$del.stable/genes.dt.all$noNa.stable
  genes.dt.all$percent.gain.o <- (genes.dt.all$gain.other-genes.dt.all$amp.other)/genes.dt.all$noNa.other
  genes.dt.all$percent.amp.o <- genes.dt.all$amp.other/genes.dt.all$noNa.other
  genes.dt.all$percent.lost.o <- (genes.dt.all$lost.other-genes.dt.all$del.other)/genes.dt.all$noNa.other
  genes.dt.all$percent.del.o <- genes.dt.all$del.other/genes.dt.all$noNa.other
  genes.dt.all[,percent.noalt.h := 1-(percent.gain.h+percent.amp.h+percent.lost.h+percent.del.h)]
  genes.dt.all[,percent.noalt.o := 1-(percent.gain.o+percent.amp.o+percent.lost.o+percent.del.o)]

  m <- match(cytoband.hg19$SYMBOL,genes.dt.all$gene)
  genes.dt.all <- genes.dt.all[na.omit(m)]
  genes.dt.all
}

percentage.all <- function(x){
  percent.loss <- (as.numeric(x[24]) + as.numeric(x[28]))/2
  percent.gain <- (as.numeric(x[22]) + as.numeric(x[26]))/2
  percent.del <- (as.numeric(x[25]) + as.numeric(x[29]))/2
  percent.amp <- (as.numeric(x[23]) + as.numeric(x[27]))/2
  data.table(percent.loss,percent.gain,percent.del,percent.amp)
}


compute.CN <- function(CN.all,annots.all,cytoband.hg19,logRatio.mean,logRatio,logRatio.TCGA,logRatio.ICGC,ref="hLMS"){
  if (ref=="hLMS"){
    CN.genes.all.noFDR <- import.copyNumber.by.genes(CN.all,annots.all,cytoband.hg19,logRatio.mean,0.01)
  }else{
    annots.oLMS <- annots.all
    annots.oLMS$Cluster <- ifelse(annots.oLMS$Cluster=="hLMS","oLMS","hLMS")
    CN.genes.all.noFDR <- import.copyNumber.by.genes(CN.all,annots.oLMS,cytoband.hg19,logRatio.mean,0.01)
  }
  genes.dt.all <- CN.genes.all.noFDR$genes.dt
  genes.dt.all <- percentage.events(genes.dt.all,cytoband.hg19)
  m <- match(genes.dt.all$gene,cytoband.hg19$SYMBOL)
  genes.dt.all$start.gene <- cytoband.hg19$gene.start[m]
  genes.dt.all <- genes.dt.all[order(start.gene,decreasing=F)]
  chrNum <- as.numeric(sub("chr","",as.vector(genes.dt.all$chr)))
  genes.dt.all <- genes.dt.all[order(chrNum,decreasing=F)]
  genes.dt.all$z.puces <- logRatio[match(genes.dt.all$gene,names(logRatio))]
  genes.dt.all$z.TCGA <- logRatio.TCGA[match(genes.dt.all$gene,names(logRatio.TCGA))]
  genes.dt.all$z.ICGC <- logRatio.ICGC[match(genes.dt.all$gene,names(logRatio.ICGC))]
  genes.dt.all$position <- 1:nrow(genes.dt.all)
  genes.dt.all$CGH <- ifelse(genes.dt.all$gene%in% row.names(CN.puces),TRUE,FALSE)
  annots.genes.all <- CN.genes.all.noFDR$annots.genes
  annots.genes.all$SYMBOL <- row.names(annots.genes.all)
  annots.genes.all <- data.table(annots.genes.all)
  colnames(annots.genes.all)[c(3,4)] <- c("z.median","fisher_genes")
  genes.dt.all$fisher_genes <- as.numeric(as.vector(annots.genes.all$fisher_genes[match(genes.dt.all$gene,annots.genes.all$SYMBOL)]))
  genes.dt.all[fisher_genes==-1 & percent.lost.h<0.25]$fisher_genes <- 0
  genes.dt.all[fisher_genes==1 & percent.gain.h<0.25]$fisher_genes <- 0
  genes.dt.all[fisher_genes==-2 & percent.del.h<0.25]$fisher_genes <- 0
  genes.dt.all[fisher_genes==2 & percent.amp.h<0.25]$fisher_genes <- 0
  genes.count.all <- genes.dt.all[,list(all=length(gene),loss=length(gene[fisher_genes==-1]),del=length(gene[fisher_genes< -1]),gain=length(gene[fisher_genes==1]),amp=length(gene[fisher_genes>1]),median.z=median(z.expr[fisher_genes!=0],na.rm=T)),by=ont]
  #genes.count.all <- genes.count.all[apply(genes.count.all[,3:6,with=F],1,function(x){sum(x!=0)})>0]
  nb.loss <- sum(genes.count.all$loss)
  nb.amp <- sum(genes.count.all$amp)
  nb.gain <- sum(genes.count.all$gain)
  nb.all <- sum(genes.count.all$all)
  nb.del <- sum(genes.count.all$del)
  fishers.band <- rbindlist(apply(genes.count.all[,2:6,with=F],1,function(x,nb.all,nb.loss,nb.amp,nb.gain,nb.del){
    f.loss <- p.adjust(fisher.test(matrix(c(x[2],x[1]-x[2],nb.loss,nb.all-nb.loss),nrow=2),alternative="greater")$p.value,n=291)
    f.amp <- p.adjust(fisher.test(matrix(c(x[5],x[1]-x[5],nb.amp,nb.all-nb.amp),nrow=2),alternative="greater")$p.value,n=291)
    f.gain <- p.adjust(fisher.test(matrix(c(x[4],x[1]-x[4],nb.gain,nb.all-nb.gain),nrow=2),alternative="greater")$p.value,n=291)
    f.del <- p.adjust(fisher.test(matrix(c(x[3],x[1]-x[3],nb.del,nb.all-nb.del),nrow=2),alternative="greater")$p.value,n=291)
    data.table(f.del,f.loss,f.gain,f.amp,enr=c(-2,-1,2,1,0)[which.min(c(f.del,f.loss,f.amp,f.gain,0.01))])
    },nb.all,nb.loss,nb.amp,nb.gain,nb.del))
  #genes.count.all <- cbind(genes.count.all,fishers.band)
  fishers.band$ont <- unique(genes.count.all$ont)
  return(list(genes.dt.all=genes.dt.all,fishers.band=fishers.band,annots.genes.all=annots.genes.all))
}


comp.rsquared <- function(compare_3studies,type){
  res <- rbindlist(apply(permutations(length(grep(type,colnames(compare_3studies))),2,grep(type,colnames(compare_3studies)),repeats.allowed=T),1,function(i,compare_3studies){
  tmp <- na.omit(as.matrix(compare_3studies[,i,with=F]))
  x <- tmp[,1]
  y <- tmp[,2]
  res <- round(summary(lm(x~y))$r.squared,2)
  res.c <- round(cor(x,y),2)
  res <- ifelse(res.c>0,res,-res)
  data.frame(g1=colnames(compare_3studies)[i[1]],g2=colnames(compare_3studies)[i[2]],RSquared=res)
  },compare_3studies))
  as.data.frame.matrix(xtabs(RSquared~g1+g2,res))
}

