format.affy.probes <- function(dir){
	affy.probes <- read.table(file.path(dir,"HG-U133_Plus_2.na34.annot.csv"),sep=",",head=TRUE)
	affy.probes <- data.table(ProbeSetID=affy.probes$Probe.Set.ID,
		gene.symbol=affy.probes$Gene.Symbol)
	affy.probes <- affy.probes[gene.symbol!="---"]
	affy.probes <- rbindlist(apply(affy.probes,1,function(x){
		genes <- unlist(strsplit(x[2]," /// "))
		tmp <- data.table(gene.symbol=genes,ProbeSetID=x[1])
		tmp <- tmp[,list(ProbeSetID,gene.symbol)]
		return(tmp)
		}))
	return(affy.probes)
}


max.probes.per.genes <- function(signal,probes){
	probes <- probes[,1:2]
	colnames(probes) <- c("ProbeSetID","gene.symbol")
	m <- match(probes$ProbeSetID,row.names(signal))
	probes <- probes[!is.na(m),]
	signal <- data.table(signal[na.omit(m),])
	signal$gene.symbol <- probes$gene.symbol
	signal <- signal[,lapply(.SD,max),by=gene.symbol]
	signal <- as.data.frame(signal)
	row.names(signal) <- signal$gene.symbol
	signal <- signal[,-1]
}

stat.probes <- function(mat1,mat2){
	genes <- row.names(mat1)
	correlation <- lapply(genes,function(gene,mat1,mat2){
			probeset1 <- mat1[gene,]
			med1 <- median(as.numeric(probeset1))
			probeset2 <- mat2[gene,]
			med2 <- median(as.numeric(probeset2))
			to.harmonize <- mean(med1,med2)
			probeset1 <- probeset1 - med1 + to.harmonize
			probeset2 <- probeset2 - med2 + to.harmonize
			probesets <- rbind(probeset1,probeset2)
			row.names(probesets) <- paste(gene,c("affy","agilent"),sep=".")
			cor.probesets <- cor(as.matrix(t(probesets)))
			all.other1 <- mat1[row.names(mat1)!=gene,]
			all.other2 <- mat2[row.names(mat2)!=gene,]
			cor.probes.set <- NULL
			for (i in 1:nrow(probesets)){
				cor.other1 <- max(cor(as.numeric(probesets[i,]),as.matrix(t(all.other1))),na.rm=TRUE)
				cor.other2 <- max(cor(as.numeric(probesets[i,]),as.matrix(t(all.other2))),na.rm=TRUE)
				cor.all <- c(cor.probesets[i,],cor.other1,cor.other2, gene)
				cor.probes.set <- rbind(cor.probes.set,cor.all)
			}	
			row.names(cor.probes.set) <- row.names(probesets)
			colnames(cor.probes.set) <- c(paste(gene,c("affy","agilent"),sep="."),"cor.other.affy","cor.other.agi","gene")
			return(cor.probes.set)

		},mat1,mat2)
		names(correlation) <- genes
		return(correlation)
}

select.genes.by.cor <- function(stats.all.genes){
	keep.genes <- unlist(lapply(stats.all.genes,function(cors){
		cors <- as.data.frame(cors)
		gene <- as.vector(unique(cors$gene))
		cors$gene <- NULL
		cors <- as.numeric(as.matrix(cors))
		cor.other <- max(cors[5:8],na.rm=T)
		same <- cors[1:2]
		cor.same <- min(same,na.rm=T)
		if (cor.same>0.8 | cor.same>cor.other){
			return(gene)		
		}
	}))
	return(keep.genes)
}

harmonize <- function(mat,n){
	affy <- mat[1:n]
	m <- n+1
	agilent <- mat[m:length(mat)]
	median.affy <- median(affy)
	median.agi <- median(agilent)
	to.harmony <- mean(c(median.affy,median.agi))
	affy <- affy - median.affy + to.harmony
	agilent <- agilent - median.agi + to.harmony
	x <- c(affy,agilent)
	return(x)
}


plot.couple.gene <- function(mat,pos1,pos2,main.title){
	mat <- as.data.frame(mat[c(pos1,pos2),])
	mat$pos <- c("pos1","pos2")
	mat <- melt(mat)
	mat$variable <- as.vector(mat$variable)
	mat$platform <- unlist(sapply(strsplit(mat$variable,"\\."),"[[",2))
	mat$sample <- unlist(sapply(strsplit(mat$variable,"\\."),"[[",1))
	mat$group <- paste(mat$pos,mat$platform,sep=".")
	g <- ggplot(mat,aes(x=sample,y=value,group=group)) + geom_line(aes(colour=group)) +
		scale_color_manual(values=c(pos1.affy="blue",pos1.agi="red",pos2.affy="green",pos2.agi="orange")) +
		theme(axis.text.x = element_blank(), panel.background = element_rect(fill = "white", colour = "black")) 
	print(g)
}


harmonize.2 <- function(gene,n){
	s <- cumsum(c(1,n[-length(n)]))
	e <- cumsum(n)
	indexes <- cbind(s,e)
	medians <- apply(indexes,1,function(x,genes){x <- median(genes[x[1]:x[2]]);return(x)},genes=gene)
	to.harmony <- median(medians)
	indexes <- cbind(indexes,medians,to.harmony)
	gene <- unlist(apply(indexes,1,function(x,genes){seg <- gene[x[1]:x[2]]; seg <- seg - x[3] + x[4] ;return(seg)},genes=gene))
	return(gene)
}

import.clinical <- function(dir,names.samp){
	#c("Type","META","differentiation","mitotic_count_10_HPF","grade_primary_tumour","sex","Localisation")
	clinical.Affy <- read.table(file.path(dir,"Annots_microarray_samples.tab"),sep="\t",head=TRUE)
	clinical.Affy[clinical.Affy==""] <- NA
	clinical.Affy <- data.table(clinical.Affy)
	#clinical.Affy[,nom_echantillon:=paste(nom_echantillon,"affy.sarcComp",sep=".")]


	#levels(clinical.Affy$Cluster) <- c("hLMS","oLMS")
	Lost_of_follow_up <- clinical.Affy$date_latest_status
	DATE=(as.numeric(as.Date(Lost_of_follow_up,format="%d/%m/%Y"))-as.numeric(as.Date(clinical.Affy$date_original_diag,format="%d/%m/%Y")))/365
	clinical.Affy[,Time_OS:=DATE]

	overallS <- ifelse(grepl("[Dd]ead",as.vector(clinical.Affy$latest_status))==T,1,0)
	clinical.Affy[,OS:=overallS]

	Lost_of_follow_up <- clinical.Affy$date_first_metastasis
	DATE=(as.numeric(as.Date(Lost_of_follow_up,format="%d/%m/%Y"))-as.numeric(as.Date(clinical.Affy$date_original_diag,format="%d/%m/%Y")))/365
	DATE[is.na(DATE)] <- clinical.Affy$Time_OS[is.na(DATE)]
	clinical.Affy[,Time_MFS:=DATE]

	clinical.Affy$Meta <- ifelse(is.na(as.vector(clinical.Affy$date_first_metastasis))==T,0,1)

	grade_primary_tumour <- clinical.Affy[,grade_primary_tumour]
	TEMP <- as.vector(grade_primary_tumour)
	TEMP[TEMP=="very low" | TEMP=="Low "] <- "1"
	TEMP[TEMP=="Intermediate "] <- "2"
	TEMP[TEMP=="High "] <- "3"
	clinical.Affy[,grade_primary_tumour:=as.numeric(TEMP)]


	clinical.Affy <- as.data.frame(clinical.Affy)
	row.names(clinical.Affy) <- clinical.Affy$nom_echantillon

	return(clinical.Affy)
}

readGSEA <- function(x,gsea.dir,res.dir,prefix){
  tables <- readHTMLTable(x,header=F)
  tables <- tables[[1]]
  if (sum(grepl("neg",x))>0){
    exps <- "neg"
  }
  else{
    exps <- "pos"
  }
  if (!is.null(tables)){
    tables <- data.table(tables[,c(2,4,6,8,11)]) 
    colnames(tables) <- c("Description","size","NES","FDR","signal")
    if (sum(grepl("GO",tables$Description))==nrow(tables)){
      exps <- paste("GO",exps,sep="_")
    } else if (sum(grepl("KEGG",tables$Description))==nrow(tables)){
      exps <- paste("KEGG",exps,sep="_")
    } else if (sum(grepl("HALLMARK",tables$Description))>0){
      exps <- paste("hallmarks",exps,sep="_")
    } else if (sum(grepl("CINSARC",tables$Description))>0){
      exps <- paste("CINSARC",exps,sep="_")
    } else if (sum(grepl("MODULE",tables$Description))>0 & sum(grepl("MORF",tables$Description))>0 & sum(grepl("GNF2",tables$Description))>0){
      exps <- paste("gene_neighbour",exps,sep="_")
    } else if (sum(grepl("UNKNOWN",tables$Description))>0){
      exps <- paste("TFBS",exps,sep="_")
    } else if (sum(grepl("REACTOME",tables$Description))>0){
      exps <- paste("datasets",exps,sep="_")
    } else if (sum(grepl("CHR[0-9]+[QP]",tables$Description))>0){
      exps <- paste("positions",exps,sep="_")
    } else {
      exps <- paste("oncogenic_signature",exps,sep="_")
    } 
    tables[,FDR:=as.numeric(as.vector(FDR))]
    tables <- tables[FDR<=0.05][order(FDR)]
    write.table(tables,file.path(res.dir,paste0(prefix,"_",exps,".tab")),sep="\t",row.names=F,quote=F)
    if(nrow(tables)>10){
    	tables <- tables[1:10]
    }
    return(tables)
  }
}

get.genes.cistarget <- function(x,sig){
		tab <- read.table(x[6],sep="\t")
		nb.genes <- length(unique(as.vector(tab[,3])))
		pc.ranked <- (nrow(tab)/tab[,1][nrow(tab)])*sig
		max.rank <- tab[,1][nrow(tab)]
		data.table(names.feature=x[1],nb.genes,pc.ranked,max.rank)
}

import.icistarget <- function(dir,group){
	up.cistarget <- readHTMLTable(file.path(dir,"report.html"),header=T,row=1)[[3]]
	names.up <- unlist(sapply(strsplit(as.vector(up.cistarget[,1]),"\n"),"[[",1))
	descr.up <- unlist(sapply(strsplit(as.vector(up.cistarget[,1]),"Description: "),"[[",2))
	descr.up <- unlist(sapply(strsplit(as.vector(descr.up),"\n"),"[[",1))
	descr.up <- unlist(sapply(strsplit(as.vector(descr.up),"[[]"),"[[",1))
	descr.up <- unlist(sapply(strsplit(as.vector(descr.up),"[(]"),"[[",1))

	cluster_PWMs_up <- read.table(file.path(dir,"clusters.tbl"),sep="\t")

	folder.up <- as.vector(up.cistarget[,7])
	fileName.up <- file.path(dir,folder.up,paste0(names.up,".targets.tsv"))
	up.cistarget <- data.table(names.feature=names.up,Description=descr.up,PWMs.cluster=cluster_PWMs_up$V1[match(names.up,cluster_PWMs_up$V3)],NES=up.cistarget$V2,type=folder.up)
	up.cistarget$fileName <- fileName.up
	up.cistarget.temp <- rbindlist(apply(up.cistarget,1,get.genes.cistarget,sig=1))
	up.cistarget.temp$group <- group
	up.cistarget <- merge(up.cistarget,up.cistarget.temp,by="names.feature")
	return(up.cistarget)
}

clinical.fisher <- function(f,df.LMS){
	fisher_clinical <- NULL
	print(f)
	feat <- df.LMS[,f]
	cl.tmp <- df.LMS$Cluster
	na <- is.na(feat)
	feat <- as.vector(feat[!na])
	cl.tmp <- as.vector(cl.tmp[!na])
	cls <- unique(cl.tmp)
	if (length(unique(feat))>2){
		for (f.tmp in unique(feat)){
			feat.temp <- feat
			feat.temp[feat.temp!=f.tmp] <- "feat2"
			feat.temp[feat.temp==f.tmp] <- "feat1"
			testfisher <- fisher.test(factor(feat.temp),factor(cl.tmp))
			tmp <- data.frame(feat1=f.tmp,feat2=paste(unique(feat[feat!=f.tmp]),collapse="+"),
				testfisher=testfisher$p.value/2,
				nb.feat1.cl1=sum(feat==f.tmp & cl.tmp==min(cls)),
				nb.feat1.cl2=sum(feat==f.tmp & cl.tmp==max(cls)),
				nb.feat2.cl1=sum(feat!=f.tmp & cl.tmp==min(cls)),
				nb.feat2.cl2=sum(feat!=f.tmp & cl.tmp==max(cls)),
				feature=f)
			fisher_clinical <- rbind(fisher_clinical,tmp)
		}
	}
	else{
		testfisher <- fisher.test(factor(feat),factor(cl.tmp))
		tmp <- data.frame(feat1=levels(factor(feat))[1],feat2=levels(factor(feat))[2],testfisher=testfisher$p.value/2,
			nb.feat1.cl1=sum(feat==levels(factor(feat))[1] & cl.tmp==min(cls)),
			nb.feat1.cl2=sum(feat==levels(factor(feat))[1] & cl.tmp==max(cls)),
			nb.feat2.cl1=sum(feat==levels(factor(feat))[2] & cl.tmp==min(cls)),
			nb.feat2.cl2=sum(feat==levels(factor(feat))[2] & cl.tmp==max(cls)),
			feature=f)
		fisher_clinical <- rbind(fisher_clinical,tmp)
	}
	fisher_clinical
}


t.score <- function(x,TCGA.annots){
	stable <- x[TCGA.annots$Cluster=="hLMS"]
	other <- x[TCGA.annots$Cluster=="oLMS"]
	score <- (mean(stable)-mean(other)) / sqrt((sd(stable)^2)/length(stable) + (sd(other)^2)/length(other))
	test <- 2*pt(-abs(score),df=(length(c(stable,other))/2)-1)
	data.table(tscore=score,ttest=test)
	}

computeCentroids <- function(annots,DF,col_interest) {
  stopifnot(length(which(colnames(annots)==col_interest))==1)
  #stopifnot(length(which(annots[,col_interest] %in% c("hLMS","oLMS")))==nrow(annots))
  stopifnot(length(which(rownames(annots) %in% names(DF)))>0)
  cl <- rep(NA,ncol(DF))
  names(cl) <- names(DF)
  nom <- names(DF)
  for (i in nom) {
    cl[i]<- annots[i,col_interest]
  }
  dd <- t(scale(t(DF),scale=F))
  L <- list()
  L$centroids <- cit.dfAggregate(dd, cl , MARGIN = 2, fAggreg = mean)
  return(L)
}

classify_LMS <- function(d2,L) {
  if (is.vector(d2)){
    N <- 1
    d2 <- na.omit(d2)
    d2 <- scale(d2)
    m <- match(row.names(L$centroids),row.names(d2))
    d2 <- cbind(d2=d2[na.omit(m)],L$centroids[!is.na(m),])
    nb.genes <- nrow(d2)
    n <- ncol(L$centroids)
    d2 <- t(d2)
    tdist <- as.matrix(cor(t(d2),method="spearman"))
    tdist <- 1-c(tdist[1,2],tdist[1,3])
    names(tdist) <- colnames(L$centroids)
    pred  <- names(tdist)[which.min(tdist)]
    tdist <- data.frame(hLMS=tdist[1],oLMS=tdist[2],pred,nb.genes=nb.genes)
  }
  else {
    N <- ncol(d2)
    d2 <- t(scale(t(d2),scale=F))
    d2=merge(d2,L$centroids,by="row.names")
    rownames(d2)=d2[,1]
    d2[,1]=NULL 
    n <- ncol(L$centroids)
    d2 <- t(d2)
    tdist <- as.matrix(cor(t(d2),method="spearman"))
    tdist <- 1-as.data.frame(tdist[1:N,(N+1):(N+n)])
    tdist$pred  <- apply(tdist,1,function(z) names(tdist)[which.min(z)])
  }
    return(tdist)
}

#cit.dfAggregate  <-  function (data, partition, MARGIN = 2, fAggreg = mean.na){
#  cMARGIN <- setdiff(c(1, 2), MARGIN)
#  n <- length(partition)
#  N <- dim(data)[MARGIN]
#  p <- dim(data)[cMARGIN]
#  if (n != N)
#    stop("Error - function cit.dfAggregate : size of partition doesn't correspond to data dimension")
#  l <- split(1:N, partition)
#  d <- data
#  if (MARGIN == 2)
#    d <- t(data)
#  d <- matrix(sapply(l, function(i) #if (length(i) == 1) {
#    unlist(d[i, ])
#  }
#  else {
#    apply(d[i, ], 2, fAggreg)
#  }), ncol = p, byrow = TRUE)
#  d <- as.data.frame(d)
#  rownames(d) <- names(l)
#  names(d) <- dimnames(data)[[cMARGIN]]
#  if (MARGIN == 2)
#    d <- as.data.frame(t(d))
#  d
#}

