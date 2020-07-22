#
# --
# JohaGL 
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(sctransform)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(cowplot)
source("~/INMG_SingleCell/scripts/functions_stock.R", local=T)
prloc = "~/INMG_SingleCell/"
prefix = "op_Filtered"  # *
resu = "results/OprescuTimePoints/"
dir.create("rds/OprescuTimePoints/", recursive = T)
daysorder = c("0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured")
dayscols = viridis(length(daysorder), alpha=0.6)
per.mito = 15
 
txtdpi.FINDER="qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_oprescuDPI.txt"
txtdpi.SCRAN="qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_oprescuDPI.txt"
txtallD0FINDER="qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_SPLITTED_4D0.txt"
txtallD0SCRAN="qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_splitted_4D0.txt"

setwd(prloc)
dir.create(resu,recursive = T)

print("* retreiving complete QualityControl info *")
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

print("intersecting doublets; copying barcodes to rownames")
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
mat = getMatrixFromtxt("data/Oprescu/oprescu_ALLraw.txt")
#correct rownames 'format' as well:
colnames(mat) = str_replace(str_replace(colnames(mat),".DPI"," DPI"),"X","")

seu <- CreateSeuratObject(mat, project="oprescuALL", min.cells = 3, min.features = 120)
head(seu@meta.data)
dim(seu@meta.data)
dim(joinedQC)

# set column $barcode :
seu@meta.data$barcode = rownames(seu@meta.data)
tmp = left_join(seu@meta.data, joinedQC, by="barcode")
rownames(tmp) = tmp$barcode
if(all(rownames(tmp) == rownames(seu@meta.data))){
  seu@meta.data <- tmp
  rm(tmp)
}else{print("something went wrong when assigning QC infos to meta.data")}

seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern="^mt-")

# ================================================================================


# ================================================================================
print("* printing quality control plots, into results *")
colstrans = c(rgb(0.54, 0.2, 0.33, 0.3), rgb(0.05, 0.27, 0.31, 0.1))
# ================================================================================
q123 <- list() ;  mycriteria=c("DOUBL_INTERSECT", "is_cell", "is_inf_outlier")
q123 <- lapply(mycriteria, function(x){
  if (x == "is_inf_outlier"){ colstrans = rev(colstrans) }
  q <- FeatureScatter(seu, "nCount_RNA", "nFeature_RNA", group.by=x, cols=colstrans) +
  scale_fill_discrete(labels = paste0(levels(seu@meta.data$x), table(seu@meta.data$x))) +
  labs(title=x)+theme(legend.text = element_text(size=8), legend.title= element_blank(), title=element_text(size=6) ) 
})

q4 <- VlnPlot(seu,features=c("nCount_RNA","nFeature_RNA","percent.mt") , split.by = "orig.ident", pt.size = 0.5)
q4g <- plot_grid(q4[[1]],q4[[2]],q4[[3]],nrow=1)
q123g <- plot_grid(plotlist = q123, nrow=1)
q1234g <- plot_grid(q123g,q4g,nrow=2)
title <- ggdraw() + draw_label(paste0("QC including scran and two doublet methods intersection"), fontface='bold')
mygeompoints = plot_grid(title, q1234g, ncol=1, rel_heights=c(0.1, 1))
pdf(paste0(resu,"QCplots.pdf"), width= 10)
mygeompoints
dev.off()

###  provisionnel: saving plots to make improvements
save(q123,colstrans, q4, file=paste0(resu,"q123colstransq4.RData"))

# ================================================================================


print("* filtering object with the criteria DOUBL_INTERSECT, is_cell, is_inf_outlier *")
# ================================================================================
filtered.seu <- subset(seu, subset= DOUBL_INTERSECT != "Doublet" &  is_cell != F &
                is_inf_outlier != T & percent.mt < per.mito)

dim(filtered.seu) #19232 51687  ==>  1506 likely doublets or NOT true cells 
saveRDS(filtered.seu)

# ================================================================================  

print("* running entire analysis, using SCTransform *")
# ================================================================================
filtered.seu = SCTransform(filtered.seu, vars.to.regress = "percent.mt", verbose = FALSE)
filtered.seu =  RunPCA(filtered.seu, verbose=FALSE)
pdf(paste0(resu,"dimHeatmap30Elbow_opres.pdf"))
DimHeatmap(filtered.seu, dims=1:30, cells = 500, balanced = TRUE)
ElbowPlot(filtered.seu, ndims = 50)
DimPlot(filtered.seu, reduction = "pca")
dev.off()

RunUMAP(filtered.seu, dims = 1:30, umap.method = "uwot")  #TODO  run fitsne as well 

filtered.seu <- FindNeighbors(filtered.seu, dims = 1:30, verbose = FALSE)
filtered.seu <- FindClusters(filtered.seu, verbose = FALSE)

filtered.seu.markersU <- FindAllMarkers(seu, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.25)

save(filtered.seu, filtered.seu.markersU, file="rds/OprescuTimePoints/filteredSEUetmarkers.RData")
# ================================================================================  

# # Interesting info: 
# # --
# # Oprescu already ejected all cells expressing less than 400 nFeatures
# # it is sad, for our next analysis it will be desirable to filter
# # using `$is_cell`and  `$is_inf_outlier` instead  :
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

