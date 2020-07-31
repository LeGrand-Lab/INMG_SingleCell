# Full Oprescu's timepoint analysis
# which incorporates MODULAR QC calculated INFO (.txt)
#     - prints QC plots into results
#     - saves filtered (and unfiltered) Seurat objects
# --
# JohaGL 
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(sctransform)
library(RColorBrewer)
library(patchwork)
library(cowplot)
source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
prloc = "~/INMG_SingleCell/"
resu = "results/OprescuTimePoints/"
rdsdir= "rds/OprescuTimePoints/"

per.mito = 15
 
txtdpi.FINDER="qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_oprescuDPI.txt"
txtdpi.SCRAN="qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_oprescuDPI.txt"
txtallD0FINDER="qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_SPLITTED_4D0.txt"
txtallD0SCRAN="qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_splitted_4D0.txt"

setwd(prloc)
dir.create(resu,recursive = T)
dir.create(rdsdir, recursive = T)
print("* retreiving complete QualityControl info (D0 and DPI)*")
# ================================================================================
print("merging D0 and DPI rows (rbind)")
allD0finder = read.table(txtallD0FINDER, sep="\t",  header=T) 
allD0scran = read.table(txtallD0SCRAN, sep="\t",header=T)
op.D0.finder = allD0finder %>% filter(sample == "OprescuD0") # pull out oprescu D0
op.D0.scran = allD0scran %>% filter(filerds == "OprescuD0") # same
rm(allD0scran,allD0finder) # remove all D0, we dont need them

op.dpi.finder =  read.table(txtdpi.FINDER, sep="\t",  header=T) 
op.dpi.scran = read.table(txtdpi.SCRAN, sep="\t",  header=T) 

op.finder = rbind(op.D0.finder, op.dpi.finder) # whole oprescu qc finder info
op.scran = rbind(op.D0.scran, op.dpi.scran) # whole oprescu qc scran info

if(print(all(as.character(op.finder$id) == as.character(op.scran$barcode)))){
  print(" all barcodes match, perfect!, setting as char type and joining")
  op.finder$id = as.character(op.finder$id)
  op.scran$barcode = as.character(op.scran$barcode)
  joinedQC = left_join(op.scran, op.finder, by=c("barcode" ="id"))
  print("if some NA introduced, the reason is: scran has ALL barcodes, whereas Finder imposed cutoff; fix the script that does DoubletsFinder step")
  # occasionnally the barcodes contain indesired "X" and "." symbols:
  joinedQC$barcode = str_replace(str_replace(joinedQC$barcode, ".DPI"," DPI"), "X","")
}else{print("oh no, barcodes do not match, stopping"); stop()}

print("intersecting doublets; passing barcodes to rownames")
dim(joinedQC); dim(joinedQC %>% distinct()) # to check dims; to check no row duplication
joinedQC <- joinedQC %>% mutate(DOUBL_INTERSECT =  case_when(
  DFclass == "Doublet" & classific == "doublets" ~ "Doublet",
  TRUE ~ "Singlet"))
# selecting only useful columns for this analysis:
joinedQC <- joinedQC %>% select(barcode, sample, DOUBL_INTERSECT, is_cell, is_inf_outlier)
rownames(joinedQC) = joinedQC$barcode
rm(op.D0.finder, op.D0.scran, op.dpi.finder, op.dpi.scran, op.finder, op.scran)
# ================================================================================  

print("* creating seurat object from raw entire matrix; QC infos to meta.data *")
# ================================================================================
mat = getMatrixFromtxt("data/Oprescu/oprescu_ALLraw.txt") # time consuming
#correct rownames 'format' as well:
colnames(mat) = str_replace(str_replace(colnames(mat),".DPI"," DPI"),"X","")
seu <- CreateSeuratObject(mat, project="oprescuALL", min.cells = 3, min.features = 120)
head(seu@meta.data) ; dim(seu@meta.data) ; dim(joinedQC)
# set column $barcode :
seu@meta.data$barcode = rownames(seu@meta.data)
tmp = left_join(seu@meta.data, joinedQC, by="barcode")  # *joinedQC is our full QCinfo*
rownames(tmp) = tmp$barcode
if(all(rownames(tmp) == rownames(seu@meta.data))){
  seu@meta.data <- tmp
  rm(tmp)
}else{print("something went wrong when assigning QC infos to meta.data")}
seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern="^mt-")
print("saving unfiltered seurat object")
saveRDS(seu, paste0(rdsdir, "NONfiltered_opre_FULL_seu.rds"))
dim(seu)#19311 53193
# ================================================================================

# ================================================================================
print("* printing quality control plots, into results *")
colstrans = c(rgb(0.54, 0.2, 0.33, 0.3), rgb(0.2, 0.66, 0.77, 0.2))
# ================================================================================
q123 <- list() ;  QCcriteria=c("DOUBL_INTERSECT", "is_cell", "is_inf_outlier")
q123 <- lapply(QCcriteria, function(x){
  tabr = table(seu@meta.data[[x]])
  captx = paste(as.character(tabr[[1]]), "/",as.character(tabr[[2]]))
  if (x != "DOUBL_INTERSECT"){ 
    colstrans = c(rgb(0.54, 0.2, 0.33, 1), rgb(0.2, 0.66, 0.77, 0.05))
        if (x == "is_inf_outlier") {colstrans = rev(colstrans) } }
  q <- FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by=x, cols=colstrans, pt.size=0.4) +
      labs(title=x, caption =captx)+
      theme(legend.text = element_text(size=9), legend.title= element_blank(), title=element_text(size=9) ) 
})
q4 <- VlnPlot(seu,features=c("nCount_RNA","nFeature_RNA","percent.mt") , split.by = "orig.ident", pt.size = 0.2)
q4g <- plot_grid(q4[[1]],q4[[2]],q4[[3]],nrow=1)
q123g <- plot_grid(plotlist = q123, nrow=1, ncol=3)
q1234g <- plot_grid(q123g,q4g,nrow=2)
title <- ggdraw() + draw_label(paste0("QC including scran and two doublet methods intersection"), fontface='bold')
mygeompoints = plot_grid(title, q1234g, ncol=1, rel_heights=c(0.1, 1))
pdf(paste0(resu,"QCplotsV.pdf"), width= 10)
mygeompoints
dev.off()
###   saving plots to make improvements , just in case
save(q123, q4, file=paste0(resu,"q123_q4_ggplots.RData"))
# ================================================================================

print("* filtering object with the criteria DOUBL_INTERSECT, is_cell, is_inf_outlier *")
# ================================================================================
filtered.seu <- subset(seu, subset= DOUBL_INTERSECT != "Doublet" &  is_cell != F &
                is_inf_outlier != T & percent.mt < per.mito)
dim(filtered.seu) #19232 51687  ==>  1506 likely doublets or NOT true cells 
saveRDS(filtered.seu, file=paste0(rdsdir,"filtered_opreFULL_seu.rds"))
rm(seu)# remove NONfiltered seu to avoid mistakes
# ================================================================================  

print("* running entire expression analysis, using SCTransform *")
# ================================================================================
filtered.seu = SCTransform(filtered.seu, vars.to.regress = "percent.mt", verbose = T)
filtered.seu =  RunPCA(filtered.seu, verbose=FALSE)
pdf(paste0(resu,"dimHeatmap30Elbow_opreFULL.pdf"))
DimHeatmap(filtered.seu, dims=1:30, cells = 500, balanced = TRUE)
ElbowPlot(filtered.seu, ndims = 50)
DimPlot(filtered.seu, reduction = "pca")
dev.off()
filtered.seu <- RunUMAP(filtered.seu, dims = 1:30)  
# using function from the file 'scripts/functions_stock.R'
filtered.seu <- KNNplusFITSNE(filtered.seu, 30, 0.8)
# FITSNE is saved into "tsne" reduction
saveRDS(filtered.seu, file=paste0(rdsdir,"filtered_opreFITSNEUMAP.rds"))
filtered.seu.markersU <- FindAllMarkers(filtered.seu, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(filtered.seu.markersU, paste0(resu,"filtered_opreMARKERSboth.txt"),
              sep="\t", col.names=T, row.names=T) 

markersTOP.pos = filtered.seu.markersU %>% filter(avg_logFC>0) %>% group_by(cluster) %>% top_n(20, wt=avg_logFC)
#as markers already ranked, no need to: '%>% arrange(p_val_adj, desc(avg_logFC), .by_group=T)' 
markersTOP.neg = filtered.seu.markersU %>% filter(avg_logFC<0) %>% group_by(cluster) %>% top_n(5, wt = -avg_logFC)
write.table(markersTOP.pos, paste0(resu,"filtered_opreMARKERS_TOPpos.txt"),
            sep="\t", col.names=T, row.names=T) 
write.table(markersTOP.neg, paste0(resu,"filtered_opreMARKERS_TOPneg.txt"),
            sep="\t", col.names=T, row.names=T) 


# ================================================================================  

# # Interesting info: 
# # --
# # Oprescu already ejected all cells expressing less than 400 nFeatures (GEO raw matrix),
# #  `$is_cell`, `$is_inf_outlier`, etc, where calculated on this splitted raw matrix . 
# # --
# test = joinedQC[joinedQC$nFeature_RNA < 400,]
# dim(test)
# [1]  0 19
# test = joinedQC[joinedQC$nFeature_RNA < 405,]
# dim(test)
# [1] 47 19
# View(test)
# table(test$is_cell)
# FALSE  TRUE 
# 2    45 
# table(test$is_inf_outlier)
# FALSE  TRUE 
# 39     8 

