suppressPackageStartupMessages(library(miRComb))
suppressPackageStartupMessages(library(citccmst))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(colorspace))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(preprocessCore))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(FactoMineR))
suppressPackageStartupMessages(library(Vennerable))
suppressPackageStartupMessages(library(mclust))
suppressPackageStartupMessages(library(XML))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(survminer))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
Txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(Txdb) <- seqlevels(Txdb)[1:24]
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(Rtsne))
suppressPackageStartupMessages(library(ggbiplot))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(RColorBrewer))


# Export current plot in different formats and sizes; from Jacques van Helden
export.plot <- function (file.prefix="PlotExport",
	                 export.formats="pdf", # supported: postscript, jpg, png, bmp, pdf
			 width=11, # in inches
			 height=8, # in inches
			 horizontal=T,
                         ... ## Additional parameters are passed to the export method
                         ) {

   ppi <- 72
   file.ext <- c(
	         postscript = "ps",
	         pdf = "pdf",
	         ps = "ps",
	         eps = "eps",
		 jpeg="jpg",
		 jpg="jpg",
		 bmp="bmp",
		 png="png",
		 svg="svg",
		 tiff="tiff")
   for (f in export.formats) {
     from.dev <- dev.cur();

     file.name <- paste(file.prefix,file.ext[f], sep=".")

     if ((f == "postscript") || (f == "ps")) {
       postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal, ...)
     } else if (f == "eps") {
       postscript(file.name,paper="special",width=width,height=height,horizontal=horizontal,onefile=F, ...)
     } else if (f == "pdf") {
       pdf(file.name, paper="special",width=width,height=height, ...)
     } else if ((f == "jpg") || (f == "jpeg")) {
       jpeg(file.name,width=(width*ppi),height=(height*ppi),quality=100, ...)
     } else if (f == "png") {
       png(file.name,width=width*ppi,height=height*ppi, ...)
     } else if (f == "bmp") {
       bitmap(file.name,width=width*ppi,height=height*ppi, ...)
     } else if (f == "svg") {
     	svg(file.name,width=width*ppi,height=height*ppi, ...)
     } else if (f == "tiff") {
     	#tiff(filename = "Rplot%03d.tiff", width = 480, height = 480, units = "px", pointsize = 12, compression = c("none", "rle", "lzw", "jpeg", "zip"), bg = "white", res = NA,  ..., type = c("cairo", "Xlib", "quartz"), antialias)
		tiff(file.name,width=width*ppi,height=height*ppi, compression = 'none', ...)
     }
      else {
       print(paste("Error: format ", f, " is not supported", sep=""))
       return()
     }
     to.dev <- dev.cur()
     dev.set(which=from.dev)
     dev.copy(which=to.dev)
     dev.set(which=to.dev)
     dev.off()
     dev.set(which=from.dev) ## This is required because dev.off() returns to the first, not the last, device
   }
}


round.pretty <- function(x, min=2){ 
  if( is.null(x) ) return(NULL)   
  d <- max(x) - min(x)
  n <- 0
  while(d<1){
    d <- d*10
    n <- n+1
  } 
  round(x, max(min,n))
}


generate.annot.col <- function(annotation){
  annots <- list()
  for (i in 1:ncol(annotation)){
    temp <- annotation[,i]
    names(temp) <- row.names(annotation)
    annots[[i]] <- temp
  }
  names(annots) <- colnames(annotation)
  annotation <- annots
  print(names(annotation))
  count = 0
  for(i in 1:length(annotation)){
    if(class(annotation[[i]]) %in% c("character", "factor")){
      # convert to factor
      if( !is.factor(annotation[[i]]) )
        annotation[[i]] <- as.factor(annotation[[i]])
      count = count + nlevels(annotation[[i]])
    }
  }
  
  palette <- c(colours()[c(26,8,142,33,139,121,57,451,226,547,578,461,24,226,53,642,6,72,474,433,419)])
    
  res_colors <- vector("list",length(annotation))
  names(res_colors) <- names(annotation)
  for(i in 1:length(annotation)){
    ann <- annotation[[i]]
    aname <- names(annotation)[i]
    acol <- NULL
    if( is.null(acol) ){
      res_colors[[i]] <-
        if(class(annotation[[i]]) %in% c("character", "factor")){
          lev <- levels(ann)
          ind = 1:length(lev)
          if (length(lev)==2){
            acol <- setNames(c(colours()[613],colours()[132]), lev)
          }
          else if (length(lev)==3){
            acol <- setNames(c(colours()[142],colours()[11],colours()[172]), lev)
          }
          else if (grepl("fisher",aname)){
              if (length(lev)==3){
                acol <- setNames(c(colours()[131],"white","red"), lev)
              }
              else {
                acol <- setNames(c(colours()[131],colours()[121],"white",colours()[373],"red"), lev)
              }
              
            }
          else{
            factor_colors <- palette[1:length(lev)]
            acol <- setNames(factor_colors, lev)
          }
          # conserve NA value               
          acol[which(is.na(names(acol)))] <- NA
          acol
        }
      else{
        h = round(runif(1) * 360)
           setNames(c(colours()[130],colours()[121],"white","red",colours()[555]), round.pretty(range(c(-2,2), na.rm=TRUE)))
      }
      
    }else{
      
      acol <- 
        if( length(acol) == 1 && grepl("^\\$", acol) ) # copy colors from other columns if the spec starts with '$'
          annotation_colors[[substr(acol, 2, nchar(acol))]]
      else if( is.factor(ann) ){
        # convert to a palette of the number of levels
        acol <- ccPalette(acol, nlevels(ann))
        
        if( is.null(names(acol)) )
          names(acol) <- levels(ann)
        
        acol
      }else{
        acol <- ccPalette(acol)
        if( is.null(names(acol)) )
          names(acol) <- round.pretty(seq(min(ann, na.rm=TRUE), max(ann, na.rm=TRUE), length.out=length(acol)))
        
        acol
      }
      
      # update the colors if necessary
      if( !is.null(acol) )
        res_colors[[i]] <- acol
    }
    
    # store type information
    attr(res_colors[[i]], 'afactor') <- is.factor(ann)    
    
  } 
  names(res_colors) <- names(annotation)
  # return ordered colors as the annotations
  return(res_colors)  
}


# From http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

import.pathways <- function(data.dir){
  hallmarks <- read.gmt(file.path(data.dir,"ext","h.all.v6.1.symbols.gmt.txt"))
  GO <- read.gmt(file.path(data.dir,"ext","c5.all.v6.1.symbols.gmt.txt"))
  positions <- read.gmt(file.path(data.dir,"ext","c1.all.v6.1.symbols.gmt.txt"))
  motifs <- read.gmt(file.path(data.dir,"ext","c3.all.v6.1.symbols.gmt.txt"))
  KEGG <- read.gmt(file.path(data.dir,"ext","KEGG.symbols.gmt.txt"))
  onc <- read.gmt(file.path(data.dir,"ext","c6.all.v6.1.symbols.gmt.txt"))
  cancer <- read.gmt(file.path(data.dir,"ext","c4.all.v6.1.symbols.gmt.txt"))
  datasets <- read.gmt(file.path(data.dir,"ext","c2.all.v6.1.symbols.gmt.txt"))
  pathways <- list(hallmarks=hallmarks,GO=GO,motifs=motifs,positions=positions,KEGG=KEGG,onc=onc,cancer=cancer,datasets=datasets)
}




