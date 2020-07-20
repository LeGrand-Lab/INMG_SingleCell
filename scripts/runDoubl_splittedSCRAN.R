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
library(DropletUtils)

prloc="~/INMG_SingleCell/"
setwd(prloc)
rdsdir = "rds/doubletsD0spli_SCRAN/"
outdir = "qcdoubl/spli_SCRAN/"
outtablename = "TABLE_DOUBLETS_SCRAN_splitted_4D0.txt"
dir.create(outdir, recursive=T)
dir.create(rdsdir, recursive=T)

#load functions : rundoublets_scran and others
source(file="~/INMG_SingleCell/scripts/functions_stock.R") 

#doSplitDeMicheli("data/DeMicheliD0/", "rawdataD0.txt") # DONE (ONLY ONCE NEEDED)

# ================================================================================
print("defining full paths raw splitted data")
# ================================================================================
dorsos <- list.files("data/DellOrsoD0", full.names=T) # subfolders inthere
gios <-  list.files("data/GiordaniD0", pattern="\\.csv$", full.names= T)
dmichs <- list.files("data/DeMicheliD0/", pattern="DeMich.txt", full.names = T)
oprescu <- list.files("data/OprescuD0", pattern="Noninjured_raw.txt$", full.names = T)
#confirmed, yes oprescu at noninjury is one single run

names(dorsos) = list.files("data/DellOrsoD0") # NOT full.names  
names(gios) <- sapply(gios, function(x) { lf = strsplit(x, "_");lf=strsplit(lf[[1]],"/")
  key = lf[[1]][ length(lf[[1]]) ] }) # like 'GSM3520458'
names(dmichs) = sapply(dmichs, function(x) { lf = strsplit(x, "/")
  key = str_replace(lf[[1]][ length(lf[[1]]) ], ".txt","")})
names(oprescu) = strsplit(oprescu,"/")[[1]][2]

ALLFILESTORUN =c(dorsos, gios, dmichs, oprescu)
typefile = c(rep("10X", length(dorsos)), rep("csv", length(gios)), rep("txt",length(dmichs)),
                    rep("txt", length(oprescu))) #  

# ================================================================================
print("Running analysis separately for each dataset")
# ================================================================================

exitstatus <- mapply(  function(x, y) {
  if (y == "10X"){
    mat = getMatrixFrom10X(ALLFILESTORUN[[x]])
  }else if (y == "csv"){
    mat = getMatrixFromCsv(ALLFILESTORUN[[x]])
  }else {
    mat =  getMatrixFromtxt(ALLFILESTORUN[[x]])
  }
  sce <- SingleCellExperiment(assays=list(counts=as.matrix(mat)))
  print(sce)
  sce <- rundoublets_scran(sce)

  print("mark empty drops")
  bcrank <- DropletUtils::barcodeRanks(SingleCellExperiment::counts(sce[rowData(sce)$expressed, ]))
  # Note: long version has $n_umi which is equal to $sum ($sum is calculated by ::nexprs)
  colData(sce)$is_cell <- colData(sce)$sum > metadata(bcrank)$inflection  # $is_cell when FALSE is an empty drop 
  print("end mark empty drops")
  
  print(paste0("saving rds into ", rdsdir, x, ".rds"))
  saveRDS(sce, paste0(rdsdir, x, ".rds"))
  
  f <- scater::plotColData(sce, x="sum",y="detected",colour_by="doublet_score")
  t <- plotTSNE(sce, colour_by="doublet_score")
  pdf(paste0(outdir, x, "plotandKnee.pdf"), width=12)
  print(plot_grid(f,t))
  print(plot_grid(knee_plot(bcrank),NULL))
  dev.off()
  return(0)
} ,names(ALLFILESTORUN), typefile)

# ================================================================================
print("creating dataframe to save all filtering info by scran, by barcode")
# ================================================================================
outdf = data.frame(barcode=character(), filerds=character(),  is_inf_outlier = logical(), is_cell = logical(), 
                   doublet_score=double(),  classific=factor())
# check rds files to run in desired order:
for(i in 1:length(ALLFILESTORUN)){
  print(file.exists(paste0(rdsdir, names(ALLFILESTORUN[i]),".rds")))
}
for(i in 1:length(ALLFILESTORUN)){
  sce <- readRDS( paste0(rdsdir, names(ALLFILESTORUN[i]),".rds") )
  dim(sce)
  tmp = data.frame(  barcode= rownames(colData(sce)),
                     filerds = rep(names(ALLFILESTORUN[i]), length(rownames(colData(sce)))),
                     is_inf_outlier = colData(sce)$keep_total,
                     is_cell = colData(sce)$is_cell, 
                     doublet_score = colData(sce)$doublet_score  )
  Q95 = quantile(sce$doublet_score,0.95)
  tmp <- tmp %>% mutate(classific = case_when(
    doublet_score >= Q95 ~ "doublets",
    TRUE ~ "singlet"  
  ))
  outdf <- rbind(outdf, tmp)
}
rm(tmp)
print(paste0("saving table (including all 4 D0 sets) into ", outdir))
write.table(outdf, paste0(outdir,outtablename), sep="\t", col.names = T)

# ================================================================================
# NOTE : example of contstruction of files paths list
# ================================================================================
# dorsos <- list.files("data/DellOrsoD0", full.names=T) # subfolders inthere
# gios <-  list.files("data/GiordaniD0", pattern="\\.csv$", full.names= T)
# dmichs <- list.files("data/DeMicheliD0/", pattern="DeMich.txt", full.names = T)
# oprescu <- list.files("data/OprescuD0", pattern="Noninjured_raw.txt$", full.names = T)
# #confirmed, yes oprescu at noninjury is one single run
# 
# names(dorsos) = list.files("data/DellOrsoD0") # NOT full.names  
# names(gios) <- sapply(gios, function(x) { lf = strsplit(x, "_");lf=strsplit(lf[[1]],"/")
# key = lf[[1]][ length(lf[[1]]) ] }) # like 'GSM3520458'
# names(dmichs) = sapply(dmichs, function(x) { lf = strsplit(x, "/")
# key = str_replace(lf[[1]][ length(lf[[1]]) ], ".txt","")})
# names(oprescu) = strsplit(oprescu,"/")[[1]][2]
# 
# ALLFILESTORUN =c(dorsos, gios, dmichs, oprescu)
# typefile = c(rep("10X", length(dorsos)), rep("csv", length(gios)), rep("txt",length(dmichs)),
#              rep("txt", length(oprescu))) #  

# ================================================================================
# ================================================================================
# OLD version with sub-classification on doublets
# :
# outdf = data.frame(barcode=character(), filerds=character(),doublet_score=double(),
#                    doubinterval=factor(), classifQ92.5=factor(), classifQ95=factor())
# # check rds files to run in desired order:
# for(i in 1:length(ALLFILESTORUN)){
#   print(file.exists(paste0(rdsdir, names(ALLFILESTORUN[i]),".rds")))
# }
# 
# for(i in 1:length(ALLFILESTORUN)){
#   sce <- readRDS( paste0(rdsdir, names(ALLFILESTORUN[i]),".rds") )
#   dim(sce)
#   tmp = data.frame(  barcode= rownames(colData(sce)),
#                      filerds = rep(names(ALLFILESTORUN[i]), length(rownames(colData(sce)))),
#                      doublet_score = colData(sce)$doublet_score  )
#   Q95 = quantile(sce$doublet_score,0.95)
#   Q92.5 = quantile(sce$doublet_score, 0.925)
#   
#   tmp <- tmp %>% mutate(classific = case_when(
#     doublet_score >= Q95 ~ "doublets high conf",
#     doublet_score > Q92.5 & doublet_score < Q95 ~ "doublets low conf",
#     TRUE ~ "singlet"  
#   ))
#   outdf <- rbind(outdf, tmp)
# }
# print(paste0("saving table into ", outdir))
# write.table(outdf, paste0(outdir,"TABLE_DOUBLETS_SCRAN_splitted.txt"), 
#             sep="\t", col.names = T)

