# Frequently used functions
# 
# --
# Joha 2020


NormFeatScalePCA <- function(seu, nFeatRNA, percentmit){
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern= "^mt-")
  seu <- subset(seu, subset = nFeature_RNA > 200 & 
                  nFeature_RNA < nFeatRNA & percent.mt < percentmit)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor=10000)
  seu <- FindVariableFeatures(seu, selection.method="vst", nfeatures= 2000)
  all.genes <- rownames(seu)
  seu <- ScaleData(seu, features = all.genes)
  seu <- RunPCA(seu,features=VariableFeatures(object = seu))  
  return(seu)
}


KNNplusFITSNE <- function(seu, nbDIM, resolu){
  #note: more resolution more groups,range 0.4-1.12
  seu <- FindNeighbors(seu, dims=1:nbDIM)
  seu <- FindClusters(seu, resolution=resolu)
  # Run FIt-SNE
  source("~/INMG_SingleCell/programsSC/FIt-SNE/fast_tsne.R", chdir=T)
  nbCELLS.div100 <- round(dim(seu)[2]/100) #nombre de cellules divisé par 100
  seu <- RunTSNE(object = seu,
                         perplexity=nbCELLS.div100, 
                         reduction="pca",
                         dims=1:nbDIM,
                         tsne.method = "FIt-SNE",
                         nthreads=4,
                         reduction.key="FItSNE_",
                         fast_tsne_path="~/INMG_SingleCell/programsSC/FIt-SNE/bin/fast_tsne",
                         max_iter=1000)
   return(seu)
}

definecolors <- function(celltype){ 
  # celltype is: seu@active.idents OR seu@meta.data$celltype
  print("for information: argument vector MUST HAVE LEVELS")
  c <- colorRampPalette(brewer.pal(8,"Set2"))(length(levels(celltype)))
  return(c)
}

customTransferLabels <- function(markersDF, refDF, NBtop)  {
  # this function takes two dataframe (markers, Areference) and an integer:
  #   (how many top positive markers from markers tibble to use)
  # return vector string having celltypes, levels corresponding to cluster numbers:
  topmarkersDF <- markersDF %>% group_by(cluster) %>% top_n(n=NBtop,wt=avg_logFC)
  joined <- left_join(topmarkersDF, refDF, by = "gene") %>%
    select("cluster","gene","clusternames","avg_logFC") %>% distinct()
  bestmatch <- joined %>% group_by(cluster) %>%
    summarize(clusternames=names(which.max(table(clusternames))))
  vec = bestmatch$clusternames
  names(vec) = bestmatch$cluster
  return(vec)
}


doCustomImputeCelltype <- function(refDF, resu, markersinresults, delimiter,
                                   rdsdir, seufile, outsuffix){
  # Function for 'impute..' scripts, SAVES celltype_Vector.rds :
  seu <-readRDS(paste0(rdsdir,seufile))
  markersDF <- read.table(paste0(resu, markersinresults ), 
                          sep=delimiter,
                          header=TRUE)
  if(length(colnames(markersDF)) < 2) {
    print("ERROR: less than two columns in markersDF, delimiter is WRONG")
    return(1)
  }else {
    matchedtypes <- customTransferLabels(markersDF,refDF, 6) #using 6top+ markers
    seu <- RenameIdents(seu, matchedtypes)
    # save 'matchedtypes', a vector compatible with ..._seu_fitsne.rds
    print("saving celltypes (vector) specific to this seurat object into rds/")
    saveRDS(matchedtypes, file=paste0(rdsdir,outsuffix))
    return(seu) # returns seurat with clusters names
  }
}


