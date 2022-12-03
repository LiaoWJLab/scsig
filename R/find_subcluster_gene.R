



#' Title
#'
#' @param sce 
#' @param assay 
#' @param col_celltype 
#' @param col_sub_celltype 
#' @param cols 
#' @param seed 
#' @param palette 
#' @param show_col 
#' @param genes_for_deg 
#' @param min_cell_count 
#' @param no_limitation_celltypes 
#' @param index 
#' @param path 
#' @param show_feas_pheatmap 
#' @param remove_other_celltypes 
#' @param find_cluster_sig 
#'
#' @return
#' @export
#'
#' @examples
find_subcluster_gene<-function(sce,
                               assay                     = "integrated",
                               col_celltype              = "celltypes",
                               col_sub_celltype          = "scpred_seurat",
                               cols                      = "normal",
                               seed                      = 123,
                               palette                   = 1,
                               show_col                  = F,
                               genes_for_deg             = NULL,
                               min_cell_count            = 30,
                               no_limitation_celltypes   = "Epithelial cells",
                               index                     = 1,
                               path                      = "4-Diff-degs-of-subclusters",
                               show_feas_pheatmap        = 8,
                               remove_other_celltypes    = TRUE,
                               find_cluster_sig          = FALSE){
  
  
  ############################################
  # col_celltype<-"celltypes"
  # col_sub_celltype<- "scpred_seurat"
  ############################################
  celltypes<-unique(sce@meta.data[,col_celltype])
  ###########################################
  
  message(">>>---Celltype of Seurat object:")
  print(table(sce[[col_celltype]]))
  
  
  for (i in index) {
    
    celltype<-celltypes[i]
    message(paste0(">>>>----Processing celltypes ", celltype))
    ############################################
    
    Idents(sce)<- col_celltype
    sces_sub<-subset(sce, idents = celltype)
    
    if(remove_other_celltypes){
      sces_sub@meta.data[, col_sub_celltype]<-ifelse(grepl(sces_sub@meta.data[, col_sub_celltype],pattern = celltype), sces_sub@meta.data[, col_sub_celltype], "unassigned" )
    }
    
    #去除细胞数量很少的clusters
    ###########################################
    if(!is.null(no_limitation_celltypes)){
      if(celltype == no_limitation_celltypes){
        ignore<- celltype
      }else{
        ignore<- NULL
      }
    }else{
      ignore<- NULL
    }
    
    sces_sub<-unassign_cell(sce            = sces_sub,
                        cluster            = col_sub_celltype,
                        ignore_cell_prefix = ignore,
                        min_cell_count     = min_cell_count,
                        new_col            = NULL,
                        delete_unassigned  = T,
                        return_meta_data   = FALSE)
    
    #' 如果目标变量没有亚组直接跳过此次循环====================
    if(length(unique(sces_sub@meta.data[, col_sub_celltype]))==1) {
      message(">>>--- cells with only one level, this celltype will be skiped...")
      print(table(sces_sub@meta.data[, col_sub_celltype]))
      next
    }
    ############################################
    
    path_res<-creat_folder(paste0(i,"-",celltype))
    ############################################
    

    res<- dong_find_markers(   sce                      = sces_sub,
                                 assay                    = assay,
                                 slot                     = "scale.data",
                                 group                    = col_sub_celltype,
                                 verbose                  = T,
                                 fig.type                 = "pdf",
                                 pt.size                  = 1,
                                 cols                     = cols,
                                 seed                     = 123,
                                 show_col                 = F,
                                 palette                  = palette,
                                 show_genes               = show_feas_pheatmap,
                                 show_genes_pheatmap      = show_feas_pheatmap,
                                 logfc.threshold          = 0.15,
                                 test.use                 = "wilcox",
                                 only.pos                 = TRUE,
                                 feature_type             = "gene",
                                 enrich_cutoff_logfc      = 0.2,
                                 enrich_cutoff_padj       = 0.05,
                                 hwidth                   = 10,
                                 hheight                  = NULL,
                                 show_features            = 8,
                                 show_plot                = T,
                                 path                     = path_res$folder_name,
                                 character_limit          = 50,
                                 recluster                = FALSE,
                                 assay_for_recluster      = NULL,
                                 dims_for_recluster       = 20,
                                 resolution_for_recluster = 0.2)
    
  }
  
  
  
  
}