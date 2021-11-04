
import.PANCAN.colors <- function(){
	X_primary_disease=c("ovarian serous cystadenocarcinoma"=colours()[8],"uterine corpus endometrioid carcinoma"=colours()[11],"uterine carcinosarcoma"=colours()[12],"cervical & endocervical cancer"=colours()[139],"pancreatic adenocarcinoma"=colours()[43],"stomach adenocarcinoma"=colours()[121],"rectum adenocarcinoma"=colours()[124],"colon adenocarcinoma"=colours()[131],"prostate adenocarcinoma"=colours()[132],"lung adenocarcinoma"=colours()[471],"lung squamous cell carcinoma"=colours()[116],"liver hepatocellular carcinoma"=colours()[369],"kidney papillary cell carcinoma"=colours()[421],"kidney clear cell carcinoma"=colours()[542],"bladder urothelial carcinoma"=colours()[463],"esophageal carcinoma"=colours()[466],"head & neck squamous cell carcinoma"=colours()[468],"cholangiocarcinoma"=colours()[455],"thyroid carcinoma"=colours()[422],"breast invasive carcinoma"=colours()[367],"glioblastoma multiforme"=colours()[228],"brain lower grade glioma"=colours()[22],"testicular germ cell tumor"=colours()[142],"sarcoma"=colours()[210],"thymoma"=colours()[656],"mesothelioma"=colours()[424],"adrenocortical cancer"=colours()[615],"kidney chromophobe"=colours()[1],"diffuse large B-cell lymphoma"=colours()[41],"pheochromocytoma & paraganglioma"=colours()[106],"skin cutaneous melanoma"=colours()[584],"uveal melanoma"=colours()[600],"hLMS"=colours()[24],"oLMS"="red")
	return(X_primary_disease)
}


run_mirComb <- function(annots,mRNA,miRNA,diff.miRNAs,DIO3){
	pheno.data <- data.frame(group=annots$Cluster,DvH=ifelse(annots$Cluster=="hLMS",1,0),row.names=row.names(annots))
	data.obj<-new("corObject",dat.miRNA=as.matrix(miRNA),dat.mRNA=as.matrix(mRNA),
                            pheno.miRNA=pheno.data,pheno.mRNA=pheno.data)
	data.obj<-addDiffexp(data.obj,"miRNA",classes="DvH",method.dif="limma",method.adj="holm")
	data.obj<-addDiffexp(data.obj,"mRNA",classes="DvH",method.dif="limma",method.adj="holm")
	data.obj<-addSig(data.obj,"mRNA",pval=0.01,logratio=1)
	data.obj<-addSig(data.obj,"miRNA",pval=0.01,logratio=1)
	data.obj<-addCorrelation(data.obj,alternative="less")
	data.obj<-addNet(data.obj)
	data.obj<-addFoldchanges(data.obj)
	data.obj<-addDatabase(data.obj,database=c("miRecords","miRTarBase"))
	data.obj<-correctPval(data.obj, pval="pval",method.adj="fdr")
	data.obj<-addScore(data.obj)
	#lapply(diff.miRNAs,function(mir,data.obj){
	#	print(mir)
	#	data.obj <<-GOanalysis(data.obj,type="REACTOME",ontology="REACTOME",pval.cutoff=0.01,dat.sum=1,sub.miRNA=mir)
	#	names(data.obj@GO.results)[names(data.obj@GO.results)=="REACTOME:REACTOME"] <- paste("REACTOME",mir,sep=":")
	#	},data.obj)
	#data.obj<-GOanalysis(data.obj,type="REACTOME",ontology="REACTOME",pval.cutoff=0.01,dat.sum=1,sub.miRNA=DIO3,up=T)
	#names(data.obj@GO.results)[names(data.obj@GO.results)=="REACTOME:REACTOME"] <- "REACTOME:DIO3"
	#data.obj<-GOanalysis(data.obj,type="REACTOME",ontology="REACTOME",pval.cutoff=0.01,dat.sum=1,up=T)
	#names(data.obj@GO.results)[names(data.obj@GO.results)=="REACTOME:REACTOME"] <- "REACTOME:UP"
	#data.obj<-GOanalysis(data.obj,type="REACTOME",ontology="REACTOME",pval.cutoff=0.01,dat.sum=1,dw=T)
	#names(data.obj@GO.results)[names(data.obj@GO.results)=="REACTOME:REACTOME"] <- "REACTOME:DN"
	# data.obj<-GOanalysis(data.obj,type="KEGG",ontology="KEGG",pval.cutoff=0.01,dat.sum=1,up=T)
	# names(data.obj@GO.results)[names(data.obj@GO.results)=="KEGG:KEGG"] <- "KEGG:KEGG:UP"
	# data.obj<-GOanalysis(data.obj,type="KEGG",ontology="KEGG",pval.cutoff=0.01,dat.sum=1,dw=T)
	# data.obj<-GOanalysis(data.obj,type="GO",ontology="BP",pval.cutoff=0.01,dat.sum=1,up=T)
	# names(data.obj@GO.results)[names(data.obj@GO.results)=="GO:BP"] <- "GO:BP:UP"
	# data.obj<-GOanalysis(data.obj,type="GO",ontology="BP",pval.cutoff=0.01,dat.sum=1,dw=T)
	# data.obj<-GOanalysis(data.obj,type="GO",ontology="CC",pval.cutoff=0.01,dat.sum=1,up=T)
	# names(data.obj@GO.results)[names(data.obj@GO.results)=="GO:CC"] <- "GO:CC:UP"
	# data.obj<-GOanalysis(data.obj,type="GO",ontology="CC",pval.cutoff=0.01,dat.sum=1,dw=T)
	# data.obj<-GOanalysis(data.obj,type="GO",ontology="MF",pval.cutoff=0.01,dat.sum=1,up=T)
	# names(data.obj@GO.results)[names(data.obj@GO.results)=="GO:MF"] <- "GO:MF:UP"
	# data.obj<-GOanalysis(data.obj,type="GO",ontology="MF",pval.cutoff=0.01,dat.sum=1,dw=T)
	return(data.obj)
}

extract.genes <- function(x){
  members <- x[8]
  members <- unlist(strsplit(members,"/"))
  members <- data.table(SYMBOL=members,term=x[2],ont=x[3],synth=x[7])
  return(members) 
}


