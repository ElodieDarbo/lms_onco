## Copy Number
## re-do ICGC
files.CN <- system(paste("ls",file.path(data.dir,"CN_by_LMS","*")),intern=T)
cn.select <- sub("/data/CN_by_LMS/CopyNumber_WL500_all_","",files.CN)
cn.select <- sub(".txt","",cn.select)
files.CN <- files.CN[cn.select%in%row.names(merge.clinical)]
cn.select <- cn.select[cn.select%in%row.names(merge.clinical)]

list.ICGC <- list()
all.genes.ICGC <- c()
for (n in 1:length(files.CN)){
	temp <- read.table(files.CN[n],sep="\t",)
	colnames(temp) <- c("seqnames","start","end","width","CN")
	temp$CN <- as.numeric(as.vector(sub("CN","",as.vector(temp$CN))))
	temp <- as(temp,"GRanges")
	ovl <- findOverlaps(genes.coords,temp)
	genes.temp <- genes.coords[queryHits(ovl)]
	genes.temp$CN <- temp$CN[subjectHits(ovl)]
	genes.temp <- data.table(as.data.frame(genes.temp,row.names=NULL))
	genes.temp <- genes.temp[,list(CN=min(CN)),by=c("seqnames","start","end","SYMBOL")]
	genes.temp <- as(as.data.frame(genes.temp),"GRanges")
	genes.temp$CN[as.vector(genes.temp$CN)>4] <- 4
	all.genes.ICGC <- unique(c(all.genes.ICGC,genes.temp$SYMBOL))
	CN.temp <- genes.temp$CN
	names(CN.temp) <- genes.temp$SYMBOL
	list.ICGC[[n]] <- CN.temp
	names(list.ICGC)[n] <- cn.select[n]
}

CN.ICGC <- as.data.frame(do.call("cbind",list.ICGC))

### TCGA
CN_TCGA <- read.table(file.path(data.dir,"TCGA_SARC_gistic2thd-2015-02-24","genomicMatrix"),sep="\t",head=T,row=1)
CN_clinical <- read.table(file.path(data.dir,"TCGA_SARC_gistic2thd-2015-02-24","clinical_data"),sep="\t",head=T,row=1)
row.names(CN_clinical) <- gsub("-","\\.",row.names(CN_clinical))
CN_TCGA <- CN_TCGA + 2

CN_clinical <- merge.clinical[row.names(merge.clinical)%in%colnames(CN_TCGA),]

CN_TCGA <- CN_TCGA[,row.names(CN_clinical)]
m <- match(row.names(CN.ICGC),row.names(CN_TCGA))
CN <- cbind(CN.ICGC[!is.na(m),],CN_TCGA[na.omit(m),])
