# DOUBLETFINDER on splitted datasets
# --
# JohaGL    

library(dplyr)
library(tidyverse)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(DoubletFinder)

prloc="~/INMG_SingleCell/"
setwd(prloc)
rdsdir = "rds/doubletsD0spli_FINDER/"
outdir = "qcdoubl/spli_FINDER/"
dir.create(outdir, recursive=T)
dir.create(rdsdir, recursive=T)
source(file="~/INMG_SingleCell/scripts/functions_stock.R", local=T)

#remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')

# ================================================================================
print("defining full paths raw splitted data, order as previous Finder run")
# ================================================================================
dorsos <- list.files("data/DellOrsoD0", full.names=T) # subfolders inthere
gios <-  list.files("data/GiordaniD0", pattern="\\.csv$", full.names= T)
dmichs <- list.files("data/DeMicheliD0/", pattern="DeMich.txt", full.names = T)
names(dorsos) = list.files("data/DellOrsoD0") # NOT full.names  
names(gios) <- sapply(gios, function(x) { lf = strsplit(x, "_");lf=strsplit(lf[[1]],"/")
    key = lf[[1]][ length(lf[[1]]) ] }) # like 'GSM3520458'
names(dmichs) = sapply(dmichs, function(x) { lf = strsplit(x, "/")
    key = str_replace(lf[[1]][ length(lf[[1]]) ], ".txt","")})
# names() : serve as output filenames !!
projectnames <- c("dorso1_wt1","dorso_wt2","gio_wt1","gio_wt2","demichili_D0_A", "demichili_D0_B", "demichili_D0_Cv3")
# project as William wrote
ALLFILESTORUN =c(dorsos,gios,dmichs)


# ================================================================================
print("running Seurat and DoubletFinder using same parameters as in prev test (test.R by William)")
# ================================================================================

metadf = data.frame(groupe=character(),nCount_RNA=integer(),nFeature_RNA=integer(),percent.mt=double(),
                    RNA_snn_res=factor(), cluster=factor(),pANN=double(),DFclass=character(),
                    id=character(),UMAP_1=double(),UMAP_2=double(),TSNE_1=double(),TSNE_2=double(),
                    sample=character())

for (i in 1:length(projectnames)){
  x = names(ALLFILESTORUN[i]); y=projectnames[i]
  print(paste0("**** ", ALLFILESTORUN[i], " **** >> out:  ",rdsdir,x ))
  if (y == "dorso1_wt1" || y == "dorso_wt2"){
    mat = getMatrixFrom10X(ALLFILESTORUN[[x]]) ; PROJ = y
  }else if (y == "gio_wt1" || y == "gio_wt2"){
    mat = getMatrixFromGio(ALLFILESTORUN[[x]]) ; PROJ = y
  }else {
    mat =  getMatrixFromtxt(ALLFILESTORUN[[x]]) ; PROJ= y 
  } 
  seu <- CreateSeuratObject(mat, project=PROJ, min.cells=3, min.features=200) 
  rm(mat)
  # **** NOTE  : trycatch concerns last sample ,
  #    D0_CvDeMich.txt, see end of code
  seu_filtered <- tryCatch({
    seu_filtered <- NormFeatScalePCA(seu, 2500, 5)
    return(seu_filtered)
  },error = function(err){
    seu_filtered <- NormFeatScalePCA(seu, 6000, 20)
  })
  rm(seu)
  seu_filtered <- FindNeighbors(seu_filtered, dims = 1:10)
  seu_filtered <- FindClusters(seu_filtered, resolution = 0.5)
  seu_filtered <- RunUMAP(seu_filtered, dims = 1:10)
  seu_filtered <- RunTSNE(seu_filtered, dims = 1:13, resolution = 0.5)
  sweep.res.seu_filtered <- paramSweep_v3(seu_filtered, PCs = 1:10, sct = FALSE)
  sweep.stats_muscle <- summarizeSweep(sweep.res.seu_filtered, GT = FALSE)
  bcmvn_muscle <- find.pK(sweep.stats_muscle)
  #table(Idents(seu_filtered), seu_filtered$clusters)
  #homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*length(seu_filtered@meta.data$orig.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  #nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seu_filtered <- doubletFinder_v3(seu_filtered, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  #seu_filtered <- doubletFinder_v3(seu_filtered, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
  
  saveRDS(seu_filtered, paste0(rdsdir,x, ".rds"))
  umapdata=data.frame(seu_filtered@reductions$umap@cell.embeddings)
  tsnedata=data.frame(seu_filtered@reductions$tsne@cell.embeddings)
  metadata=seu_filtered@meta.data
  umapdata$id=as.character(rownames(umapdata))
  metadata$id=as.character(rownames(metadata))
  tsnedata$id=as.character(rownames(tsnedata))
  plotdoublet = inner_join(x=metadata,y=umapdata,by="id")
  plotdoublet = inner_join(x=plotdoublet,y=tsnedata,by="id")
  rm(umapdata,metadata)
  colnames(plotdoublet) = c("groupe","nCount_RNA","nFeature_RNA","percent.mt",
                            "RNA_snn_res","cluster","pANN","DFclass","id",
                            "UMAP_1","UMAP_2","TSNE_1","TSNE_2")
  plotdoublet$DFclass[plotdoublet$DFclass %in% "Doublet" & plotdoublet$pANN <= 0.5]="Doublets - Low confidence"
  plotdoublet$DFclass[plotdoublet$DFclass %in% "Doublet" & plotdoublet$pANN > 0.5 ]="Doublets - High confidence"
  
  # new col: sample=> needed for prefixes:
  plotdoublet$sample = rep( y,  length(plotdoublet$id) )
  metadf = rbind(metadf, plotdoublet)
} 

write.table(metadf, paste0(outdir, "TABLE_DOUBLFINDER_SPLITTED.txt"), sep="\t", col.names=T)






# ================================================================================

# IMPORTANT NOTE   : concerning 'D0_CvDeMich' sample 
# *** DEBUGGING INITIALLY WAS **:
# seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern= "^mt-")
# seu <- subset(seu, subset = nFeature_RNA > 200 &
#                 nFeature_RNA < 2500 & percent.mt < 5)
# seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor=10000)
# seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures= 2000)
# all.genes <- rownames(seu)
# seu <- ScaleData(seu, features = all.genes)
# seu_filtered <- RunPCA(seu,features=VariableFeatures(object = seu), npcs=5) # ==> tested several,npcs=5 finally ok

# tested also CSTransform, did not work, in conclusion, THE PROBLEM WAS THE STRINGENT FILTERS, 
# fixed as authors: up to 6000 feat admitted, up to 20% mitochondrial  :


# for (i in 7:length(projectnames)){
#   print(paste0("**** ", ALLFILESTORUN[i], " **** >> out:  ",rdsdir,x ))
#   x = names(ALLFILESTORUN[i]); y=projectnames[i]
#   if (y == "demichili"){
#     mat =  getMatrixFromtxt(ALLFILESTORUN[[x]]) ; PROJ = y
#   }else if (y == "gio_wt1" || y == "gio_wt2"){
#     mat = getMatrixFromGio(ALLFILESTORUN[[x]]) ; PROJ = y
#   }else {
#     mat = getMatrixFrom10X(ALLFILESTORUN[[x]]) ; PROJ = y
#   } 
#   seu <- CreateSeuratObject(mat, project=PROJ, min.cells=3, min.features=200) 
#   
#   seu_filtered <- tryCatch({
#     seu_filtered <- NormFeatScalePCA(seu, 2500, 5)
#     print("partie try")
#   },error = function(err){
#     seu_filtered <- NormFeatScalePCA(seu, 6000, 20)
#   })
#   rm(seu)
#   seu_filtered <- FindNeighbors(seu_filtered, dims = 1:10)
#   seu_filtered <- FindClusters(seu_filtered, resolution = 0.5)
#   seu_filtered <- RunUMAP(seu_filtered, dims = 1:10)
#   seu_filtered <- RunTSNE(seu_filtered, dims = 1:13, resolution = 0.5)
#   sweep.res.seu_filtered <- paramSweep_v3(seu_filtered, PCs = 1:10, sct = FALSE)
#   sweep.stats_muscle <- summarizeSweep(sweep.res.seu_filtered, GT = FALSE)
#   bcmvn_muscle <- find.pK(sweep.stats_muscle)
#   #table(Idents(seu_filtered), seu_filtered$clusters)
#   #homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
#   nExp_poi <- round(0.075*length(seu_filtered@meta.data$orig.ident))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
#   #nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
#   seu_filtered <- doubletFinder_v3(seu_filtered, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
#   #seu_filtered <- doubletFinder_v3(seu_filtered, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
#   
#   saveRDS(seu_filtered, paste0(rdsdir,x, ".rds")) # 
#   
#   umapdata=data.frame(seu_filtered@reductions$umap@cell.embeddings)
#   tsnedata=data.frame(seu_filtered@reductions$tsne@cell.embeddings)
#   metadata=seu_filtered@meta.data
#   umapdata$id=as.character(rownames(umapdata))
#   metadata$id=as.character(rownames(metadata))
#   tsnedata$id=as.character(rownames(tsnedata))
#   plotdoublet = inner_join(x=metadata,y=umapdata,by="id")
#   plotdoublet = inner_join(x=plotdoublet,y=tsnedata,by="id")
#   rm(umapdata,metadata)
#   colnames(plotdoublet) = c("groupe","nCount_RNA","nFeature_RNA","percent.mt",
#                             "RNA_snn_res","cluster","pANN","DFclass","id",
#                             "UMAP_1","UMAP_2","TSNE_1","TSNE_2")
#   plotdoublet$DFclass[plotdoublet$DFclass %in% "Doublet" & plotdoublet$pANN <= 0.5]="Doublets - Low confidence"
#   plotdoublet$DFclass[plotdoublet$DFclass %in% "Doublet" & plotdoublet$pANN > 0.5 ]="Doublets - High confidence"
#   metadf = rbind(metadf, plotdoublet)
# }

# ================================================================================



