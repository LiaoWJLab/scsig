



#'  Find and show gene expression markers for cell subtypes
#'
#'  Find and show markers (differentially expressed genes) for cell subytpes in a dataset
#'
#' @param sce Seurat object with cell type annotation.
#' @param assay Assay to pull from, e.g. RNA, SCT, integrated
#' @param col_celltype Name of column where cell types in the metadata
#' @param col_sub_celltype Name of column where cell subtypes in the metadata
#' @param cols Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param show_col Whether to show color palettes
#' @param min_cell_count Minimal cell counts in each cluster
#' @param no_limitation_celltypes A character vector corresponding to cell types with no cell count limitation
#' @param index User can choose specific cell types to save computing source, default is 1
#' @param path Path of the output saving directory
#' @param show_feas_pheatmap Number of genes displayed in the heatmap. Default is 8
#' @param remove_other_celltypes Default is TRUE
#' @param feas Genes to test. Default is to use all genes
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.15 Increasing logfc.threshold speeds up the function, but can miss weaker signals
#' @param fig.type Format of plot saving, such as pdf and png
#' @param pt.size Size of point
#' @param show_feas_dim Number of genes displayed in the violin plot, default is 8
#' @param recluster If TRUE, assay_for_recluster must be provided
#' @param assay_for_recluster Name of assay to use, such as RNA, SCT, integrated, default is integrated
#' @param re_scale_tsne_umap
#'
#' @return Seurat object containing specific type of cells
#' @export
#'
#' @examples
find_subcluster_gene<-function(sce,
                               assay                     = "integrated",
                               slot                      = "scale.data",
                               col_celltype              = "celltypes",
                               col_sub_celltype          = "scpred_seurat",
                               feas                      = NULL,
                               cols                      = "normal",
                               seed                      = 123,
                               palette                   = 1,
                               show_col                  = F,
                               min_cell_count            = 30,
                               no_limitation_celltypes   = "Epithelial cells",
                               index                     = 1,
                               path                      = "4-Diff-degs-of-subclusters",
                               show_feas_pheatmap        = 8,
                               show_feas_dim             = 12,
                               remove_other_celltypes    = TRUE,
                               logfc.threshold           = 0.15,
                               fig.type                  = "pdf",
                               pt.size                   = 0.7,
                               re_scale_tsne_umap        = FALSE,
                               recluster                 = F,
                               assay_for_recluster       = "integrated"){


  ############################################
  # col_celltype<-"celltypes"
  # col_sub_celltype<- "scpred_seurat"
  ############################################
  celltypes<-unique(sce@meta.data[,col_celltype])
  ############################################

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
    }else{
      sces_sub@meta.data[, col_sub_celltype]<-ifelse(grepl(sces_sub@meta.data[, col_sub_celltype],pattern = celltype), sces_sub@meta.data[, col_sub_celltype], "other_celltypes" )
    }

    #Remove the clusters with few cells
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

    #' If the target variable has no subgroup, skip the loop directly====================
    if(length(unique(sces_sub@meta.data[, col_sub_celltype]))==1) {
      message(">>>--- cells with only one level, this celltype will be skiped...")
      print(table(sces_sub@meta.data[, col_sub_celltype]))
      next
    }
    ############################################

    path_res<-creat_folder(paste0(path,"/", i,"-",celltype))
    ############################################


         dong_find_markers(    sce                       = sces_sub,
                                assay                    = assay,
                                slot                     = slot,
                                feas                     = feas,
                                group                    = col_sub_celltype,
                                verbose                  = T,
                                fig.type                 = fig.type,
                                pt.size                  = pt.size,
                                cols                     = cols,
                                seed                     = 123,
                                show_col                 = F,
                                palette                  = palette,
                                show_genes               = show_feas_pheatmap,
                                show_genes_pheatmap      = show_feas_pheatmap,
                                logfc.threshold          = logfc.threshold,
                                test.use                 = "wilcox",
                                only.pos                 = TRUE,
                                feature_type             = "gene",
                                enrich_cutoff_logfc      = 0.2,
                                enrich_cutoff_padj       = 0.05,
                                hwidth                   = 10,
                                hheight                  = NULL,
                                show_features            = show_feas_dim,
                                show_plot                = T,
                                path                     = path_res$folder_name,
                                character_limit          = 50,
                                re_scale_tsne_umap        = re_scale_tsne_umap,
                                recluster                = recluster,
                                assay_for_recluster      = assay_for_recluster,
                                dims_for_recluster       = 20,
                                resolution_for_recluster = 0.2)

  }

  return(sces_sub)

}
