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

prloc = "~/INMG_SingleCell/"
prefix = "op_Filtered"  # *
resu = "results/OprescuTimePoints"
daysorder = c("0.5 DPI", "2 DPI", "3.5 DPI", "5 DPI", "10 DPI", "21 DPI", "Noninjured")
dayscols = viridis(length(daysorder), alpha=0.6)
 
txtdpi.FINDER="qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_oprescuDPI.txt"
txtdpi.SCRAN="qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_oprescuDPI.txt"
txtallD0FINDER="qcdoubl/spli_FINDER/TABLE_DOUBLFINDER_SPLITTED_4D0.txt"
txtallD0SCRAN="qcdoubl/spli_SCRAN/TABLE_DOUBLETS_SCRAN_splitted_4D0.txt"

setwd(prloc)
dir.create(resu,recursive = T)

print("retreiving complete QualityControl info, merging D0 and DPI rows (rbind)")
# ================================================================================
allD0finder = read.table(txtallD0FINDER, sep="\t",  header=T) 
allD0scran = read.table(txtallD0SCRAN, sep="\t",header=T)
op.D0.finder = allD0finder %>% filter(sample == "OprescuD0")
op.D0.scran = allD0scran %>% filter(filerds == "OprescuD0")
rm(allD0scran,allD0finder)

op.dpi.finder =  read.table(txtdpi.FINDER, sep="\t",  header=T) 
op.dpi.scran = read.table(txtdpi.SCRAN, sep="\t",  header=T) 

op.finder = rbind(op.D0.finder, op.dpi.finder)
op.scran = rbind(op.D0.scran, op.dpi.scran)

if(print(all(as.character(op.finder$id) == as.character(op.scran$barcode)))){
  print(" all barcodes match, perfect!, setting as char type")
  op.finder$id = as.character(op.finder$id)
  op.scran$barcode = as.character(op.scran$barcode)
}else{print("oh no, barcodes do not match, stopping"); stop()}

joinedQC = left_join(op.finder, op.scran, by=c("id" = "barcode"))
dim(joinedQC); dim(joinedQC %>% distinct())

joinedQC <- joinedQC %>% mutate(DOUBL_INTERSECT =  case_when(
  DFclass == "Doublet" & classific == "doublets" ~ "Doublet",
  TRUE ~ "Singlet"))
joinedQC <- joinedQC %>% select(id, sample, DOUBL_INTERSECT, is_cell, is_inf_outlier)
#### CORRECT BARCODES!!!!!  replace .DPI. by " DPI_" and X by "" 
rownames(joinedQC) = 
  temp = str_replace(joinedQC$id
print("creating seurat object from raw entire matrix, add QC to metadata")
# ================================================================================

# ================================================================================  
# ================================================================================

print("printing quality control plots, into results")
# ================================================================================

# ================================================================================



print("")
# ================================================================================

# ================================================================================  

print("")
# ================================================================================

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

