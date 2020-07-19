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
rdsdir = "rds/doubletsD0spli_SCRAN/"
outdir = "qcdoubl/spli_SCRAN/"
dir.create(outdir, recursive=T)
dir.create(rdsdir, recursive=T)
#load functions : rundoublets_scran and others
source(file="~/INMG_SingleCell/scripts/functions_stock.R") 

#doSplitDeMicheli("data/DeMicheliD0/", "rawdataD0.txt") # DONE (ONLY ONCE NEEDED)

# ================================================================================
print("defining full paths raw splitted data")
# ================================================================================
dmichs <- list.files("data/DeMicheliD0/", pattern="DeMich.txt", full.names = T)
gios <-  list.files("data/GiordaniD0", pattern="\\.csv$", full.names= T)
dorsos <- list.files("data/DellOrsoD0", full.names=T) # subfolders inthere
names(dmichs) = sapply(dmichs, function(x) { lf = strsplit(x, "/")
    key = str_replace(lf[[1]][ length(lf[[1]]) ], ".txt","")})
names(gios) <- sapply(gios, function(x) { lf = strsplit(x, "_");lf=strsplit(lf[[1]],"/")
            key = lf[[1]][ length(lf[[1]]) ] }) # like 'GSM3520458'
names(dorsos) = list.files("data/DellOrsoD0") # NOT full.names  

ALLFILESTORUN =c(dorsos, gios, dmichs)


# ================================================================================
print("Running analysis separately for each dataset")
# ================================================================================
exitstatus <- lapply( names(ALLFILESTORUN), function(x) {
    if (x == "dorsowt1" | x == "dorsowt2"){
      mat = getMatrixFrom10X(ALLFILESTORUN[[x]])
     }else if (x == "GSM3520458" | x == "GSM3520459"){
      mat = getMatrixFromGio(ALLFILESTORUN[[x]])
    }else {
      mat =  getMatrixFromtxt(ALLFILESTORUN[[x]])
    }
    sce <- SingleCellExperiment(assays=list(counts=as.matrix(mat)))
    print(sce)
    sce <- rundoublets_scran(sce)
    scores = as.numeric(as.character(sce$doublet_score))
    print(paste0("saving rds into ", rdsdir, x, ".rds"))
    saveRDS(sce, paste0(rdsdir, x, ".rds"))
    f <- scater::plotColData(sce, x="sum",y="detected",colour_by="doublet_score")
    t <- plotTSNE(sce, colour_by="doublet_score")
    pdf(paste0(rdsdir, x, "plot.pdf"), width=12)
    print(plot_grid(f,t))
    dev.off()
  return(0)
  } )

exitstatus <- lapply( names(ALLFILESTORUN), function(x) {
  if (x == "dorsowt1" | x == "dorsowt2"){
    #mat = getMatrixFrom10X(ALLFILESTORUN[[x]])
  }else if (x == "GSM3520458" | x == "GSM3520459"){
    #mat = getMatrixFromGio(ALLFILESTORUN[[x]])
  }else {
    #mat =  getMatrixFromtxt(ALLFILESTORUN[[x]])
  }
  print(x)
  #sce <- readRDS(paste0("saving rds into ", rdsdir, x, ".rds"))
  
  return(0)
} )


# ================================================================================
print("creating dataframe : $barcode $filerds $doubletscore $classific")
# ================================================================================
# then open -One by one- the  sce objects :

outdf = data.frame(barcode=character(), filerds=character(),doublet_score=double(),
                   doubinterval=factor(), classifQ92.5=factor(), classifQ95=factor())
# check rds files to run in desired order:
for(i in 1:length(ALLFILESTORUN)){
  print(file.exists(paste0(rdsdir, names(ALLFILESTORUN[i]),".rds")))
}

for(i in 1:length(ALLFILESTORUN)){
  sce <- readRDS( paste0(rdsdir, names(ALLFILESTORUN[i]),".rds") )
  dim(sce)
  tmp = data.frame(  barcode= rownames(colData(sce)),
                   filerds = rep(names(ALLFILESTORUN[i]), length(rownames(colData(sce)))),
                doublet_score = colData(sce)$doublet_score  )
  Q95 = quantile(sce$doublet_score,0.95)
  Q92.5 = quantile(sce$doublet_score, 0.925)

  tmp <- tmp %>% mutate(classific = case_when(
    doublet_score >= Q95 ~ "doublets high conf",
    doublet_score > Q92.5 & doublet_score < Q95 ~ "doublets low conf",
    TRUE ~ "singlet"  
  ))

  outdf <- rbind(outdf, tmp)

}

print(paste0("saving table into ", outdir))
write.table(outdf, paste0(outdir,"TABLE_DOUBLETS_SCRAN_splitted.txt"), 
            sep="\t", col.names = T)

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

