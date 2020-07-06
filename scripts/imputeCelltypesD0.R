#!/usr/bin/env Rscript
# Impute cell types to D0 independent datasets
# input: 
#   txt file ==> all specific markers from experience
#   txt file => reference markers
#   rds file 
# output : .rds file containing corresponding celltypes (vector)
#         and figures
# --
# JohaGL
library(dplyr)
library(Seurat)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
source(file="~/INMG_SingleCell/scripts/functions_stock.R",local=T)

prloc="~/INMG_SingleCell/"
setwd(prloc)

#          *-* Vars *-*
delimiter = " " # espaced delimited results/../ALLMARKERS_....txt
#   Markers groups and celltypes as discussed with Dr Le Grand:
ref = "refmarkers/newRefmarkersToCells_v1.txt" # >>check
refDF <- read.table(ref,sep="\t",header=TRUE)
head(refDF,n=2)
#clusternames   gene
#1 B_lymphocytes.CD4   Cd74
#2 B_lymphocytes.CD4  Cd79a
#      

outsuffix="celltype_Vector.rds"
delimiter = " " # space delimited results/../ALLMARKERS_....txt
# define an order for running the loop:
authors = c("DellOrsoD0", "DeMicheliD0",
              "GiordaniD0","OprescuD0")
# define vector .rds files SAME ORDER and assign names:
rds.name = c("dorso_seu_fitsne.rds","dmizero_seu_fitsne.rds",
             "gio_seu_fitsne.rds","opre_seu_fitsne.rds")
names(rds.name) = authors

# ================================================================================
# imputing cell types to each D0 :
# ================================================================================

print("checking vector as you defined it, to run analysis")
print(rds.name)

seulist <- list()
seulist <- lapply(names(rds.name), function(p) doCustomImputeCelltype(
                              refDF,
                              paste0("results/",p,"/"),
                              paste0("ALLMARKERS_",p,".txt"),
                              delimiter, 
                              paste0("rds/",p,"/"),
                              as.character(rds.name[p]),
                              outsuffix)
         )
print("generate DimPlot of each seurat named object")
print("using mapply to parse 'seulist' and 'rds.name' matching indexes") 
printorder <- mapply(function(seu,p){
  print(seu)
  pdf(paste0("results/",p,"/", "FITSNE_named.pdf"))
  plo <- DimPlot(seu, label=T,repel=T, 
          cols=definecolors(seu@active.ident))+
    ggtitle(str_replace(p,"D0"," D0"))
  print(plo) # MUST BE to print into file as needed 
  dev.off()
  return(p)}, seulist, names(rds.name)
)
print("order in which plots were saved :")
print(printorder)
print("")
print("END")
print("In subsequent jobs, to get seurat object with its respective cell type, 
      open both (seu.rds and vector.rds) and use Seurat function 'RenameIdents'")

# END
# ================================================================================

#  Previous version not using loop (dangerous, prone to errors):
# # ================================================================================
# # imputing cell types to Oprescu D0:
# # ================================================================================
# resu4 = "results/OprescuD0/"  # >> check
# rdsdir4 = "rds/OprescuD0/" # >> check
# seufile4 = "opre_seu_fitsne.rds"   # >> check
# markersinresults4 = "ALLMARKERS_OprescuD0.txt"   # >> check
# 
# oprescutypeseu <- doCustomImputeCelltype(refDF, resu4, markersinresults4, delimiter,
#                                        rdsdir4, seufile4, outsuffix)
# 
# pdf(paste0(resu4, "FITSNE_named.pdf"))
# DimPlot(oprescutypeseu, label=T,repel=T, 
#         cols=definecolors(oprescutypeseu@active.ident))+
#   ggtitle("DellOrso D0")
# dev.off()

# ================================================================================
# === NOTE :
# ================================================================================
# 'doCustomImputeCelltype'  does this (no function version):

# markersDF <- read.table(paste0(resu, markersinresults ), 
#                         sep=delimiter,
#                         header=TRUE)
# seu <-readRDS(paste0(rdsdir,seufile))
# 
# matchedtypes <- customTransferLabels(markersDF,refDF, 6) #using 6top+ markers
# print(tail(matchedtypes,n=3))
# #  15                        16                        17 
# # "B_lymphocytes.CD4" "Neutrophyls_macrophages"               "Myonuclei" 
# 
# # check cell types have been imputed already or not
# unique(seu@active.ident)
# # impute cell types using 'matchedtypes' vector
# seu <- RenameIdents(seu, matchedtypes)
# unique(seu@active.ident)  # we see only 13 levels because homonyms are merged
# 
# mycols <- colorRampPalette(brewer.pal(8,"Set2"))(length(levels(seu@active.ident)))
# 
# pdf(paste0(resu, "FITSNE_named.pdf"))
# DimPlot(seu, label=T,repel=T, cols=mycols)
# dev.off()
# 
# # save 'matchedtypes', a vector compatible with ..._seu_fitsne.rds
# saveRDS(matchedtypes, file=paste0(rdsdir,outsuffix))

