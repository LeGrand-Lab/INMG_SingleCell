# Run doubletCells (scran) for each experiment separately
# ALSO: split DeMicheliD0 into separated experiments !!
# note that Giordani and DellOrso exist in splitted GEO elements.
# saved into data.
# --
# JohaGL
library(dplyr)
library(tidyverse)
library(Seurat)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(RColorBrewer)
library(cowplot)


prloc="~/INMG_SingleCell/"
setwd(prloc)
outdir = "rds/doubletsD0spli_sce/"
dir.create(outdir, recursive=T)

# ================================================================================

doSplitDeMicheli <- function(datadirD, filenm){
  dmizerol <- read.table(paste0(datadirD, filenm), sep="\t",
                         header=T, row.names=1)
  dim(dmizerol)
  # see experiments on barcodes (colnames) prefixes:
  myexps = unique(substring(colnames(dmizerol), first=1, last=5)) #"D0_A_" "D0_B_" "D0_Cv"
  listofexps = sapply(myexps, function(x) {
    select(dmizerol, starts_with(x))
  })
  listofexps <- setNames(listofexps,myexps)
  nbcolssum = 0
  for (i in 1:length(listofexps)){ nbcolssum = nbcolssum + dim(listofexps[[i]])[2]}
  
  if( (dim(dmizerol)[2]) == nbcolssum){
    print(paste0("split operation sucessful, saving into ", datadirD))
    lapply( names(listofexps), function(k) { 
      write.table( listofexps[[k]],
        paste0(datadirD, k,"DeMich.txt") , sep="\t", col.names = T);return(0)
    })
  }else{ print("errors in dimensions of splitted tables, check, nothing to save")
    return(1)}
}

getMatrixFromtxt <- function(filepath){
  m = read.table(filepath, sep="\t", header=T, row.names = 1)
  return(as.matrix(m))
}
getMatrixFrom10X <- function(dirpath){
  m = Read10X(data.dir=dirpath)
  return(as.matrix(m))
}
getMatrixFromGio <- function(filepath){
  m = read.csv( filepath, sep=",", header=TRUE, row.names=1)
  return(as.matrix(m))
}

rundoublets <- function(sce){
  rowData(sce) <- DataFrame(genes_names = rownames(sce))
  rowData(sce)$expressed <- scater::nexprs(sce,byrow=TRUE)>0
  per_cell <- perCellQCMetrics(sce[rowData(sce)$expressed,])
  colData(sce) <- cbind(colData(sce),per_cell)
  head(rowData(sce))
  colData(sce)$keep_total <- scater::isOutlier(colData(sce)$sum,type = "lower", log=TRUE)
  table(colData(sce)$keep_total) # TRUE are OUTLIERS
  sce <- scater::addPerFeatureQC(sce)
  print("Finding doublets")
  sce <- computeSumFactors(sce) # by scran: scaling normalization (implements deconvolution)
  sce <- logNormCounts(sce)
  dbl_dens <- doubletCells(sce)
  sce$doublet_score <- 0
  sce$doublet_score <- log10(dbl_dens + 1)
  sce <- runTSNE(sce, perplexity=50, PCA=T, num_threads=4)
  return(sce)
  }

# ================================================================================
### run:

#doSplitDeMicheli("data/DeMicheliD0/", "rawdataD0.txt")

# define full paths :
dmichs <- list.files("data/DeMicheliD0/", pattern="DeMich.txt", full.names = T)
gios <-  list.files("data/GiordaniD0", pattern="\\.csv$", full.names= T)
dorsos <- list.files("data/DellOrsoD0", full.names=T) # subfolders inthere
names(dmichs) = sapply(dmichs, function(x) { lf = strsplit(x, "/")
    key = str_replace(lf[[1]][ length(lf[[1]]) ], ".txt","")})
names(gios) <- sapply(gios, function(x) { lf = strsplit(x, "_");lf=strsplit(lf[[1]],"/")
            key = lf[[1]][ length(lf[[1]]) ] }) # like 'GSM3520458'
names(dorsos) = list.files("data/DellOrsoD0") # NOT full.names  

ALLFILESTORUN =c(dmichs,gios,dorsos)
prin("defining type of dataset (author), same order as in ALLFILESTORUN")
typesauth <- c(rep("dmich",length(dmichs)) ,  rep("gio",length(gios)) ,
               rep("dorso",length(dorsos)) )


print("Running analysis separately for each dataset")
x <- mapply( function(x,y){
    if (y == "dmich"){
     mat =  getMatrixFromtxt(ALLFILESTORUN[[x]])
    }else if (y == "gio"){
      mat = getMatrixFromGio(ALLFILESTORUN[[x]])
    }else {
      mat = getMatrixFrom10X(ALLFILESTORUN[[x]])
    }
    sce <- SingleCellExperiment(assays=list(counts=as.matrix(mat)))
    print(sce)
    sce <- rundoublets(sce)
    scores = as.numeric(as.character(sce$doublet_score))
    saveRDS(sce, paste0(outdir, x, ".rds"))
    f <- scater::plotColData(sce, x="sum",y="detected",colour_by="doublet_score")
    t <- plotTSNE(sce, colour_by="doublet_score")
    pdf(paste0(outdir, x, "plot.pdf"), width=12)
    print(plot_grid(f,t))
    dev.off()
  return(0)
  }, names(ALLFILESTORUN), typesauth )

# then open -One by one- the seven sce objects and do intervals
# to create a dataframe : $barcode $author $doubletscore $doubinterval
myrds = list.files("rds/doubletsD0spli_sce",pattern="\\.rds$", full.names = T)

outdf = data.frame(barcode=character(), author=character(),doublet_score=double(),doubinterval=factor(),
                   classificat=factor())

for(i in 1:length(myrds)){
  sce <- readRDS(myrds[i])
  dim(sce)
  tmp = data.frame(  barcode= rownames(colData(sce)),
                   author = rep(typesauth[i], length(rownames(colData(sce)))),
                doublet_score = colData(sce)$doublet_score  )
  
  top = quantile(sce$doublet_score, 0.995)
  mid = quantile(sce$doublet_score,0.95)
  
  tmp <- tmp %>% mutate(doubinterval = case_when(
                          doublet_score >= top ~ paste0(as.character(round(top,1)), " to ", 
                                                    as.character(round(max(doublet_score)),1), " :q ",names(top)),
                          doublet_score > mid & doublet_score < top ~ paste0(
                            as.character(round(mid,1))," to ", as.character(round(top,1)) ," :q ",names(mid)),
                          TRUE ~ paste0(round(mid,1),"  and less") 
                        ))
  
  tmp <- tmp %>% mutate(classific = case_when(
    doublet_score >= top ~ "doublets high conf",
    doublet_score > mid & doublet_score < top ~ "doublets low conf",
    TRUE ~ "singlet"  
  ))
  
  colData(sce)$doubinterval = tmp$doubinterval
  pdf(paste0(outdir, i, "_TEST.pdf"))
  plotTSNE(sce,colour_by = "doubinterval")
  dev.off()
  plotReducedDim(sce, dimred="TSNE", colour_by = "doubinterval")
  outdf <- rbind(outdf, tmp)
  outdf$doubinterval = as.factor(outdf$doubinterval)
}

write.table(outdf, paste0(outdir,"TABLE_DOUBLETS_splitted.txt"), col.names = T)

# ==========older version consumes too much memory:
# # list of the different experiments
# getMatricesFromtxt <- function(filepath){
#   lfilepath <- as.list(filepath)
#   # define items names : suffixes of each file, like "D0_CvDeMich"
#   names(lfilepath) <- lapply(lfilepath, function(x) { 
#     lf = strsplit(x, "/")
#     key = str_replace(lf[[1]][ length(lf[[1]]) ], ".txt","")})
#   #print(lfilepath)
#   listofmatrices <- lapply(lfilepath, function(x){
#     #m = read.table(x, sep="\t", header=T, row.names = 1)
#     df = array(c(3,6,4,7), dim=c(2,2))
#     colnames(df) = c("a","b")
#     rownames(df) = c("x", "y")
#     return(df)
#   }) 
#   return(listofmatrices)
# }
# getMatricesFrom10X <- function(dirpath){
#   print("function getMatricesFrom10X expects subfolders named differently")
#   subfo = as.list( list.files(dirpath, full.names = T) )
#   names(subfo) = list.files(dirpath) # NOT full.names
#   listofmatrices = lapply(subfo, function(x){
#     df = array(c(3,6,4,7), dim=c(2,2))
#     colnames(df) = c("a","b")
#     rownames(df) = c("x", "y")
#     return(df)
#   })
#   return(listofmatrices)
# }
# getMatricesFromGio <- function(filepath){
#   lfilepath <- as.list(filepath)
#   names(lfilepath) <- lapply(lfilepath, function(x) { 
#     lf = strsplit(x, "_"); key = lf[[1]][1] }) # like 'GSM3520458'
#   gios <- lapply(lfilepath, function(x){
#     df = array(c(3,6,4,7), dim=c(2,2))
#     colnames(df) = c("a","b")
#     rownames(df) = c("x", "y")
#     #read.csv( x, sep=",", header=TRUE, row.names=1)
#     return(df)
#   })
#   return(gios)
# }

# dmichs <- getMatricesFromtxt(c("data/DeMicheliD0/D0_A_DeMich.txt", "data/DeMicheliD0/D0_B_DeMich.txt",
#                                "data/DeMicheliD0/D0_CvDeMich.txt"))
# gios <- getMatricesFromGio(c("GSM3520458_20171018_uninjured_wt_filtered.csv",
#                              "GSM3520459_20180917_uninjured_wt_filtered.csv"))
# dorsos <- getMatricesFrom10X("data/DellOrsoD0")

