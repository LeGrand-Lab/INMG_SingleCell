# Run doubletCells (scran) for each experiment separately
# ALSO: split DeMicheliD0 into separated experiments !!
# note that Giordani and DellOrso exist in splitted GEO elements.
# saved into data.
# --
# JohaGL
library(dplyr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)


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

# list of the different experiments
getMatricesFromtxt <- function(vecfilepath){
  dmichs <- list( sapply(vecfilepath, function(x){
    print(x)
    lf = base::strsplit(x, "/")
    print(lf[3])
    key = str_replace(lf[length(lf)], ".txt","")
    print(key)
    #m = read.table(x, sep="\t", header=T, row.names = 1)
    key = c(6,7)
  }) )
  return(dmichs)
}
getMatricesFromCsv <- function(vecfilepath){
  read.csv( mmmm, 
           sep=",", header=TRUE, row.names=1)
}
getMatricesFrom10X <- function(dirpath){
  ###PENDING
}

rundoublets(){}
# ================================================================================

### run:

doSplitDeMicheli("data/DeMicheliD0/", "rawdataD0.txt")

dmichs <- getMatricesFromtxt(c("data/DeMicheliD0/D0_A_DeMich.txt", "data/DeMicheliD0/D0_B_DeMich.txt",
                                "data/DeMicheliD0/D0_CvDeMich.txt"))
gios <- getMatricesFromCsv(c("GSM3520458_20171018_uninjured_wt_filtered.csv",
            "GSM3520459_20180917_uninjured_wt_filtered.csv"))
dorsos <- getMatricesFrom10X()

SingleCellExperiment(assays=list(counts=matdat))


