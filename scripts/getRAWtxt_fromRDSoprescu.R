# Get raw .txt tables from Oprescu
# input: .rds file, which is an 'AnnotatedDataFrame' type object
# output: 
#     - oprescu_Noninjured_raw.txt  
# if desired, activate other tables (DPI or ALL) txt
library("dplyr")
library(Seurat)

prloc="~/INMG_SingleCell/" # << check
KEYNAME="oprescu" # <<< check
datafile="GSE138826_regen_data.rds"  # <<< check

setwd(prloc)

dir.create("data/Oprescu",recursive = T)
print("to run this , download oprescu raw rds file: 'GSE138826_regen_data.rds' 
      from GEO into folder 'data/Oprescu'")
KEYNAME <- gsub("/",'',KEYNAME)
datafile <- gsub("/",'', datafile)
oprefull <- readRDS(paste0("data/Oprescu/",datafile))
# exploring:
dim(as.data.frame(oprefull[["RNA"]]@data))  
# [1] 19311 53193
head(rownames(as.data.frame(oprefull[["RNA"]]@data)))
# [1] "Xkr4"    "Sox17"   "Mrpl15"  "Lypla1"  "Gm37988" "Tcea1" 

getNI <- function(oprefull){
  noninjured <- as.data.frame(oprefull[["RNA"]]@data) %>% select(starts_with("Noninjured"))
  dim(noninjured)
  dir.ni = "data/OprescuD0"
  # [1] 19311  5670
  print(paste0("saving non injured into ", dir.ni))
  write.table(noninjured, paste0(dir.ni, "/",KEYNAME,"_Noninjured_raw.txt"),
              sep="\t", col.names=T,row.names = T)
  return(0)
}

getDPI <- function(oprefull){
  dpi <- as.data.frame(oprefull[["RNA"]]@data) %>% select(!starts_with("Noninjured"))
  write.table(dpi, paste0("data/Oprescu/",KEYNAME, "_DPI_raw.txt"), 
              sep="\t", col.names=T,row.names = T)
  return(0)
}

getALL <- function(oprefull){
  write.table(as.data.frame(oprefull[["RNA"]]@data), 
              paste0("data/Oprescu/", KEYNAME, "_ALLraw.txt"),
              sep="\t", col.names=T, row.names = T)
  return(0)
}

#get only non injured
#getNI(oprefull)
getALL(oprefull)
# activate to get the other table(s):
# getDPI(oprefull)

