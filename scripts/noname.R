#!/usr/bin/env Rscript 
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
library(DoubletFinder)

prloc="~/INMG_SingleCell/"
setwd(prloc)
rdsdir = "rds/doubletsD0spli_SCRAN/"
outdir = "qcdoubl/spli_SCRAN/"
dir.create(outdir, recursive=T)
dir.create(rdsdir, recursive=T)
#load functions : rundoublets_scran and others
source(file="~/INMG_SingleCell/scripts/functions_stock.R") 


rm(ALLFILESTORUN) # make sure no previous vector not expected
oprescu = list.files("data/OprescuD0", pattern="Noninjured_raw.txt$", full.names = T)
#confirmed, yes oprescu at noninjury is one single run
names(oprescu) = strsplit(oprescu,"/")[[1]][2]
ALLFILESTORUN = c(oprescu)
typefile = c(rep("txt", length(oprescu))) #  ,rep("10X", length(dorsos)), rep("csv", length(gios) ....(dmichs).)

sce <- readRDS("rds/doubletsD0spli_SCRAN/OprescuD0.rds")
dim(sce)#[1] 19311  5670

# ================================================================================
print(" dataframe PLUSOPRESCU : $barcode $filerds $doubletscore $classific")
print(" the output has new classification based on 5% expected doublets")
# ================================================================================
# reopen saved dataframe
outdf <- read.table(paste0(outdir,"TABLE_DOUBLETS_SCRAN_splitted.txt"), 
                    sep="\t", header = T)
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
  tmp <- tmp %>% mutate(classific = case_when(
    doublet_score >= Q95 ~ "doublet",
    TRUE ~ "singlet"  
  ))
  outdf <- rbind(outdf, tmp)
}
print(paste0("saving table (including all 4 D0 sets) into ", outdir))
write.table(outdf, paste0(outdir,"TABLE_DOUBLETS_SCRAN_splitted_4D0.txt"), sep="\t", col.names = T)

# ================================================================================
#  DOUBLET FINDERRRRR
# ================================================================================
print("running DoubletFinder")
rm(outdir, rdsdir) # redefine rdsdir and outdir
rdsdir = "rds/doubletsD0spli_FINDER/"
outdir = "qcdoubl/spli_FINDER/"
metadf = read.table( paste0(rdsdir, "TABLE_DOUBLFINDER_SPLITTED.txt"), sep="\t", header=T)

for (i in 1:length(typefile)){
  x = names(ALLFILESTORUN[i]); y = typefile[i]
  print(paste0("**** ", ALLFILESTORUN[i], " **** >> out:  ",rdsdir,x ))
  if (y == "10X"){
    mat = getMatrixFrom10X(ALLFILESTORUN[[x]]) ; PROJ = x
  }else if (y == "csv"){
    mat = getMatrixFromCsv(ALLFILESTORUN[[x]]) ; PROJ = x
  }else {
    mat =  getMatrixFromtxt(ALLFILESTORUN[[x]]) ; PROJ= x 
  } 
  seu <- CreateSeuratObject(mat, project=PROJ, min.cells=3, min.features=200) 
  rm(mat)
  # **** NOTE  : including more features both Oprescu and DeMicheli
  seu_filtered <- tryCatch({
    seu_filtered <- NormFeatScalePCA(seu, 6000, 15)
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
  nExp_poi <- round(0.05*length(seu_filtered@meta.data$orig.ident))  ## CHANGED: (5%)
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
  # new col: sample=> needed for prefixes:
  plotdoublet$sample = rep( x,  length(plotdoublet$id) )
  metadf = rbind(metadf, plotdoublet)
} 

# PROVISIONNEL: MIX DOUBLETS AS "doublets"
metadf <- metadf %>% mutate(DFclass = case_when(DFclass != "Singlet" ~ "Doublet", TRUE ~ "Singlet"))
write.table(metadf, paste0(outdir, "TABLE_DOUBLFINDER_SPLITTED_4D0.txt"), sep="\t", col.names=T)



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(VennDiagram)
library(dplyr)
library(RColorBrewer)

scranres <- read.table("qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_splitted_4D0.txt",
                       sep="\t", header=T)
finderres <- read.table("qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_SPLITTED_4D0.txt", 
                        sep="\t", header=T)
print("it is normal if scran results has many more cells (no seurat filters applied yet)")
paste("scran= ", dim(scranres)[1],"  Finder= ",dim(finderres)[1])

# if problems , fix to obtain homologated prefixes to barcodes:
#ATTENTION: barcodes on demicheli have already D0_A_ , or _D0_B_ ,  or D0_Cv3_

scranres <- scranres %>% mutate(newid = case_when(
  filerds == "dorsowt1" ~ 
    paste0(as.character(filerds),"_",as.character(barcode)),
  filerds == "dorsowt2" ~ 
    paste0(as.character(filerds),"_",as.character(barcode)),
  TRUE ~ as.character(barcode)
))
finderres <- finderres %>% mutate(newid = case_when(
  groupe == "dorso1_wt1" ~ paste0("dorsowt1","_",as.character(id)),
  groupe == "dorso_wt2" ~ paste0("dorsowt2","_",as.character(id)),
  TRUE ~ as.character(id)
))
print("checking that no duplicated rows exist ( distinct() ):")
dim(finderres %>% distinct()); dim(scranres %>% distinct())
comparedf <- left_join(finderres, scranres, by="newid")

print("equivalences in classification columns")
unique(comparedf$classific)
unique(comparedf$DFclass)

vennliB <- list(
  "scran_DOUBLET" = comparedf$newid[comparedf$classific == "doublets"] ,
  "Finder_DOUBLET" = comparedf$newid[comparedf$DFclass == "Doublet"] ,
  "scran_singlet" = comparedf$newid[comparedf$classific == "singlet"] ,
  "Finder_singlet" = comparedf$newid[comparedf$DFclass == "Singlet"]
)


b <- venn.diagram(vennliB, fill=2:5, alpha=0.3, filename=NULL, margin=0.1)
grid.newpage()

print("getting barcodes in the intersection")
trueids_doubl = intersect(vennliB[["scran_DOUBLET"]], vennliB[["Finder_DOUBLET"]])

LIKELYDOUBL <- comparedf %>% dplyr::filter(newid %in% trueids_doubl)

c <- ggplot(LIKELYDOUBL, aes(x=filerds, y=pANN)) + geom_boxplot() + 
  labs(title=paste0("DOUBLETS (intersection Finder and scran methods) n= ", dim(LIKELYDOUBL)[1]))
c2 <- ggplot(LIKELYDOUBL, aes(x=filerds, y=doublet_score)) + geom_boxplot()
print("doublets proportion in each sample:")
prop = table(LIKELYDOUBL$filerds) / ( table(comparedf$filerds) / 100 )
tprop = dim(LIKELYDOUBL)[1]/(dim(comparedf)[1]/100)
# D0_A_DeMich D0_B_DeMich D0_CvDeMich    dorsowt1    dorsowt2  GSM3520458  GSM3520459 
# 0.2702703   1.9438445   3.8647343   0.7522837   1.3066202   2.4089936   2.2005295 
d <- ggplot(as.data.frame(prop), aes(Var1,Freq)) + geom_col(fill="turquoise4") +  
  labs(y="proportion of calculated doublets by sample (%)", caption = paste("total %=",tprop)) + coord_flip()

umapinter1 <- ggplot(LIKELYDOUBL, aes(x=UMAP_1,y=UMAP_2, label=cluster)) + geom_point(aes(colour=factor(classific))) + 
  labs(title="scran::doubletCells") + theme(legend.position="none") + geom_text(check_overlap=T, size=3, alpha=0.7)
umapinter2 <- ggplot(LIKELYDOUBL, aes(UMAP_1,UMAP_2, label=cluster)) + geom_point(aes(colour=factor(DFclass))) +
  labs(title="DoubletFinder") + theme(legend.title=element_blank()) + geom_text(check_overlap=T, size=3, alpha=0.7)

clus1 <- ggplot(comparedf, aes(x=cluster,fill=classific)) + geom_bar(position="fill")+
  theme(legend.position="none")
clus2 <- ggplot(comparedf, aes(x=cluster,fill=DFclass)) + geom_bar(position="fill")


pdf("qcdoubl/FINDERintersectSCRAN_4D0.pdf", width=10)
plot_grid(b, ncol=1)
plot_grid(c,c2,d, ncol = 1, rel_widths = c(4,4,2))
plot_grid(umapinter1,umapinter2, clus1 , clus2, ncol=2, rel_widths = c(3,4))
dev.off()



