# Get raw .txt tables from Oprescu
# input: .rds file, which is an 'AnnotatedDataFrame' type object
# output: 
#     - oprescu_Noninjured_raw.txt   - oprescu_ALL.txt
#     - extracted Meta Data   - splitted tables:
#               "...0.5.DPI" "...2DPI_" ... "...21DPI"
library(tidyverse)
library(Seurat)

prloc="~/INMG_SingleCell/" # << check
KEYNAME="oprescu" # <<< check
datafile="GSE138826_regen_data.rds"  # <<< check
datadir = "data/Oprescu/"
KEYNAME <- gsub("/",'',KEYNAME)
datafile <- gsub("/",'', datafile)

setwd(prloc)
dir.create(datadir,recursive = T)

print("to run this , download oprescu raw rds file: 'GSE138826_regen_data.rds' 
      from GEO into folder 'data/Oprescu'")

oprefull <- readRDS(paste0(datadir,datafile))
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

getALL <- function(oprefull){
  write.table(as.data.frame(oprefull[["RNA"]]@data), 
              paste0(datadir, KEYNAME, "_ALLraw.txt"),
              sep="\t", col.names=T, row.names = T)
  return(0)
}

getMetadata <- function(oprefull){
  md <- oprefull[["RNA"]]@varMetadata
  write.table(md, paste0(datadir,KEYNAME,"_ExtractedMetaData.txt"),
              sep="\t", col.names=T,row.names = T) 
  return(0)
}

splitDPI <- function(tabfull){
  # fix colnames
  oldcols = colnames(tabfull)
  colnames(tabfull) = str_replace(str_replace(oldcols,".DPI", " DPI"),"X","")
  print(unique(substring(colnames(tabfull),1,7)))
  patterns = c("Noninjured", "0.5 DPI_", "2 DPI_", "3.5 DPI_", "5 DPI_", "10 DPI_", "21 DPI_")
  # getDPI(tableoprefull)
  patternsDPI = patterns[patterns != "Noninjured"]
  
  sapply(patternsDPI, function(x){
    print(x)
    tmp <- tabfull %>% select(starts_with(x))
    nfile = str_replace(x," DPI_","DPI")
    write.table(tmp, paste0(datadir,"oprescu_",nfile, ".txt"), sep="\t", row.names=T, col.names=T)
  })
  return(0)
}

#getNI(oprefull) #get only non injured
#getALL(oprefull)
getMetadata(oprefull)
tabfull = read.table(paste0(datadir, KEYNAME, "_ALLraw.txt"),sep="\t",  row.names = 1,  header=T)
splitDPI(tabfull)

