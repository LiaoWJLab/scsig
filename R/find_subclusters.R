




#' Title
#'
#' @param sce
#' @param celltype
#' @param combined_subclusters
#' @param graph.name
#' @param resolution
#' @param delete_new_variable
#' @param filtering_mini_cell
#' @param min_cell_count
#' @param ignore_cell_prefix
#'
#' @return
#' @export
#'
#' @examples
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
  #' 批量寻找个细胞大类的亚组细胞=================================

  for(i in 1:length(celltypes)){
    celltype<-celltypes[i]
    sce<-FindSubCluster(sce,
                        cluster = celltype,
                        graph.name = graph.name,
                        subcluster.name = subcluster_names[i],   #增加一个新的列到数据中
                        resolution = resolution,
                        algorithm = 2)
  }
  #########################################

  #'将以上细胞亚组整合到一个列中======combined_column==================================
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
    #'去掉一些离散的细胞
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



