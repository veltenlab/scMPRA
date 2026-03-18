require(Seurat)
require(ggplot2)
require(tidyverse)
require(scmap)
require(SingleCellExperiment)
require(parallel)
require(plyr)

#Script to project on any Seurat reference, example: Triana et al., Nature Immunology 2021.
#Let s be a seurat object you wish to project.
#healthy_reference <- readRDS(url("https://figshare.com/ndownloader/files/28408638")) #better download data before...
#s <- project_anyref(s, ref = healthy_reference, ref_reduction = "MOFAUMAP", ref_pseudotime = "Myelocytes")
#script will add a new dimensionality reduction to s, called projected, in the coordinates of ref_reduction.
#It will also add three metadata columns to s:
#projected.cluster : Projected Idents() of the reference
#projected.score : Similarity score
#projected.pseudo : Projected pseudotime, where the reference pseudotime is contained in the metadata field passed as ref_pseudotime

#this is what i want to run 
#which genes are used here? maybe i also want to specify my TAP genes here... 
test <- project_anyref(chelseas, neighbours = 5, cores = 6, features = blub,
                       integrated, ref_reduction = "umapscan", ref_pseudotime = NULL,
                       save.maps = NULL, save.final = NULL)

seurat = chelseas
neighbours = 5
cores = 6
features = blub 
ref = integrated_trp53
ref_reduction = "umap"
ref_pseudotime = NULL 
save.maps = NULL
save.final = NULL 


project_anyref <- function(seurat, neighbours = 5, cores = 6, features = NULL,
                           ref, ref_reduction = "MOFAUMAP", ref_pseudotime = "Myelocytes",
                           save.maps = NULL, save.final = NULL){
  
  # Normalise data before projecting
  # seurat <- seurat %>% NormalizeData()
  
  # get normalised counts
  norm_counts <- seurat@assays$RNA@data #this was with the old seurat, with umapscan this works also 
  #norm_counts <- seurat@assays$SCT@data #this is with seurat5
  
  # load seurat object from CloneTracer manuscript
  
  
  if(is.null(features)){
    
    features <- VariableFeatures(ref)#!!
    
  }
  
  # only use variable features from AML map, otherwise the object is too large
  data_projection <- norm_counts[intersect(rownames(norm_counts), features), ]
  
  # find nearest neighbour in the map 
  sce_Culture <- SingleCellExperiment(assays = list(normcounts =  as.matrix(data_projection)))
  logcounts(sce_Culture) <- normcounts(sce_Culture)
  rowData(sce_Culture)$feature_symbol <- rownames(sce_Culture)
  
  #subset 
  ref = subset(ref, features = features)
  ref_mat = as.matrix(ref@assays$RNA@layers$data)
  rownames(ref_mat) = rownames(ref)
    
  #sce_All <- SingleCellExperiment(assays = list(normcounts = as.matrix(ref@assays$RNA@layers$data[features,]))) idk this didnt work 
  sce_All <- SingleCellExperiment(assays = list(normcounts = ref_mat))
  colnames(sce_All) <- colnames(ref)
  
  logcounts(sce_All) <- normcounts(sce_All)
  # use gene names as feature symbols
  rowData(sce_All)$feature_symbol <- rownames(sce_All)
  # remove features with duplicated names
  sce_All <- sce_All[!duplicated(rownames(sce_All)), ]
  
  done <- F  
  if(!is.null(save.maps)) {
    if(file.exists(save.maps)) {
      done <-T
      Culture_Map <- readRDS(save.maps)
    }
  }
  
  if (!done) {
    sce_Culture<-setFeatures(sce_Culture,features =  rownames(sce_Culture))
    sce_Culture <- indexCell(sce_Culture)
    
    sce_All<-setFeatures(sce_All,features =  rownames(sce_All))
    sce_All <- indexCell(sce_All)
    
    Culture_Map <- scmapCell(
      projection = sce_Culture,
      index_list = list(
        sce_All = metadata(sce_All)$scmap_cell_index
      ),
      w = neighbours)
  }
  
  if(!is.null(save.maps)) {
    if(!file.exists(save.maps)) {
      saveRDS(Culture_Map, save.maps)
    }
  }
  
  
  #calculate Coordinates
  #id = colnames(Culture_Map$sce_All[[1]])[1] #for debugging 
  #cult = Culture_Map$sce_All[[1]] #for debugging
  
  Calc<-function(id,cult){
    u <- cult[,id]
    xcoords <- ref@reductions[[ref_reduction]]@cell.embeddings[,1][u]
    ycoords <- ref@reductions[[ref_reduction]]@cell.embeddings[,2][u]
    x=median(xcoords)
    y=median(ycoords)
    x_mean=mean(xcoords)
    y_mean=mean(ycoords)
    meandist <- mean(sqrt((xcoords -x)^2 + (ycoords - y)^2))
    sdx <- sd(xcoords-x)
    sdy <- sd(ycoords-y)
    ct.t=table(Idents(ref)[u])
    ct.t<- ct.t[order(ct.t,decreasing = T)]
    Prop<-prop.table(ct.t)
    Prop<- Prop[order(Prop,decreasing = T)]
    subset_nearest <- subset(ref, features = features, cells = u) #i added this im not sure 
    nearest <- subset_nearest@assays$RNA@layers$data
    query <- normcounts(sce_Culture)[,id]
    #anglenn <- apply(nearest,2,function(x) angle(x, query))
    cornn<- apply(nearest,2,function(x) cor(x, query))
    if (is.null(ref_pseudotime)) {
      pst.projected <- NA
    } else {
      pst <- seurat@meta.data[u,ref_pseudotime]
      pst.projected<-mean(na.omit(pst))  
    }
    
    data.frame(row.names = id, x = x, y =y,x_mean=x_mean,y_mean=y_mean, ct = names(ct.t)[1],prop=Prop[1],meandist=meandist, sdx=sdx, sdy=sdy, cor = mean(cornn),pseudo_myel = pst.projected)
  }
  
  
  done <- F  
  if(!is.null(save.final)) {
    if(file.exists(save.final)) {
      done <-T
      mapped <- readRDS(save.final)
    }
  }
  
  if (!done) {
    mapped <- mclapply(colnames(Culture_Map$sce_All[[1]]), Calc, cult = Culture_Map$sce_All[[1]], mc.cores = cores) #here i get error, amazing! and how do i fix it? 
    mapped <- do.call(rbind,mapped) %>% rownames_to_column(var = "cell_barcode") %>% 
      dplyr::select(cell_barcode, x, y, ct, cor, pseudo_myel) %>% 
      dplyr::mutate(type = "query")
    
    colnames(mapped) <- c("cell_barcode", "umapx", "umapy", "projected.cluster", "projected.score", "projected.pseudo", "projected.type")
  }
  
  
  if(!is.null(save.final)) {
    if(!file.exists(save.final)) {
      saveRDS(mapped, save.final)
    }
  }
  
  
  embeddings <-as.matrix(mapped[, c("umapx","umapy")])
  rownames(embeddings) <- mapped$cell_barcode
  colnames(embeddings) <- c("PROJECTED_1","PROJECTED_2")
  seurat@reductions$projected <- CreateDimReducObject(embeddings = embeddings, key = "PROJECTED_", assay = "RNA")
  
  newmeta <- mapped[,c("projected.cluster","projected.score","projected.pseudo")]
  rownames(newmeta) <- mapped$cell_barcode
  seurat <- AddMetaData(seurat, newmeta)
  
  
  return(seurat)
}


ref_pseudotime = NULL




