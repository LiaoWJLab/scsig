




#' Sub-cluster determination
#' 
#' Identify sub-clusters for each cluster of cells
#' @param sce Seurat object
#' @param celltype Name of the metadata column to group cells by
#' @param combined_subclusters New name of a metadata column, dedault is "celltype_subcluster"
#' @param graph.name Name of graph to use for the clustering algorithm
#' @param resolution Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities
#' @param delete_new_variable Dedault is TRUE
#' @param filtering_mini_cell Whether to filter out the clusters with few cells, dedault is TRUE
#' @param min_cell_count Minimum number of cell counts of each cluster, dedault is 50
#' @param ignore_cell_prefix A character vector corresponding to cell types with no cell count limitation
#'
#' @return Suerat object with updated metadata containing subgroups of cells
#' @export
#'
#' @examples
#' ## Not RUN
#' data("pbmc_small")
#' sce<-find_subclusters(sce=pbmc_small, celltype= "RNA_snn_res.1")
#' ## Not RUN
find_subclusters<-function(sce, celltype, combined_subclusters = "celltype_subcluster",
                           graph.name = NULL, resolution = 0.3, delete_new_variable = T,
                           filtering_mini_cell = TRUE, min_cell_count = 50,ignore_cell_prefix = "null"){
  
  
  message(paste0("+++++++++++++++++++++++++++++++++++++++++++++++"))
  message(paste0(">>>== Enshure that running following code before this process....."))
  message(paste0('sce <- FindNeighbors(object = sce, assay = "SCT", reduction = "harmony", dims = 1:50, graph.name = "sct")'))
  message('sce <- FindClusters(object = sce, resolution = 0.6)')
  message('sce <- RunTSNE(sce, reduction = "harmony", dims = 1:50) ')
  message(paste0("+++++++++++++++++++++++++++++++++++++++++++++++"))
  
  Idents(sce)<- celltype
  message(">>>== Celltype Propotion: ")
  print(table(as.character(sce@meta.data[,celltype])))
  #####################################
  sce@meta.data[,combined_subclusters] <- sce@meta.data[, celltype]
  #####################################
  celltypes<-unique(Idents(sce))
  subcluster_names<- gsub(celltypes, pattern = " ", replacement = "_")
  subcluster_names<- gsub(subcluster_names, pattern = "-", replacement = "_")
  #######################################
  
  message(paste0("Graph name that can be choosed:  ", names(sce@graphs)))
  #' Batch search for subgroups of cells of a group of cell =================================
  
  for(i in 1:length(celltypes)){
    celltype<-celltypes[i]
    sce<-FindSubCluster(sce,
                        cluster = celltype,
                        graph.name = graph.name,
                        subcluster.name = subcluster_names[i],   
                        resolution = resolution,
                        algorithm = 2)
  }
  #########################################
  
  #'Integrating the above cell subgroups information into one column======combined_column==================================
  for (i in 1:length(celltypes)) {
    var<-celltypes[i]
    cellname<-gsub(var, pattern = " ", replacement = "_")
    cellname<-gsub(cellname, pattern = "-", replacement = "_")
    
    subcluster_name<-subcluster_names[i]
    sce@meta.data[grep(sce@meta.data[,combined_subclusters], pattern = var), combined_subclusters]<-
      sce@meta.data[grep(sce@meta.data[,combined_subclusters], pattern = var), subcluster_name]
    
    if(delete_new_variable){
      sce@meta.data<-sce@meta.data[,-which(colnames(sce@meta.data)==cellname)]
    }
  }
  
  print(head(sce))
  if(filtering_mini_cell){
    
    message(paste0("+++++++++++++++++++++++++++++++++++++++++++++++"))
    message(paste0(">>>== Filtering celltyes....."))
    #####################################
    #'Removal of some discrete cells
    sce<-unassign_cell(sce =  sce,
                       cluster = combined_subclusters,
                       ignore_cell_prefix = ignore_cell_prefix,
                       min_cell_count = min_cell_count,
                       new_col = NULL,
                       delete_unassigned = T,
                       return_meta_data = FALSE)
  }
  
  return(sce)
}



