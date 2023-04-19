


#' Cell type annotation by scType
#'
#' Cell type annotation using specific marker combinations from single-cell transcriptomic data
#' refer to [ScType](https://www.nature.com/articles/s41467-022-28803-w)
#' @param sce Seurat object
#' @param assay Assay to pull from, e.g. RNA, SCT, integrated
#' @param slot Data slot to use, choose from 'counts', 'data', or 'scale.data'
#' @param scale Whether the matrix is scaled，default is NULL
#' @param cluster A vector of variables to group cells by
#' @param point.size Size of point, default is 1.5
#' @param db_ Database of manually collected cell type annotation, default is "ScTypeDB_full.xlsx" deposited in data
#' @param cols Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param show_col Whether to show color palettes
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param gs_sctyper Default is null, user can provide gene set file manually
#' @param tissue_type Default is NULL
#' @param cell_type Cell types choosed from "base", "epithelial", "myeloid", "tcell", "bcell", "fibroblast" and "endothelial", default is "base"
#' @param cell_subset Default is NULL
#'
#' @return Seurat object with updated metadata containing cell type annotation
#' @export
#'
#' @examples
#' tnbc<-sctyper_anno(sce = sce, tissue = "Immune system")
sctyper_anno<-function(sce,
                       gs_sctyper = NULL,
                       tissue_type= NULL,
                       cell_type  = "base",
                       cell_subset= NULL,
                       study      = NULL,
                       assay      = NULL,
                       slot       = "scale.data",
                       scale      = NULL,
                       cluster    = "seurat_clusters",
                       db_        = "ScTyperDB-merged.xlsx",
                       db_path    = NULL,
                       gene_names_to_uppercase = TRUE){


  if(!is.null(db_path)){
    db_<- db_path
  }else{
    db_<- paste0(base::system.file("data", package = "scsig"),"/", db_)
  }

  # prepare gene sets
  gs_sctyper_list = gene_sets_prepare(gs_sctyper = gs_sctyper, path_to_db_file = db_, cell_type = cell_type, tissue_type = tissue_type, cell_subset = cell_subset)
  ###################################

  message(">>>---Assay used to estimation:")
  if(!is.null(assay)){
    print(DefaultAssay(sce))
  }else{
    print(paste0(">>>>> ",assay))
  }


  # get cell-type by cell matrix
  scRNAseqData<- SeuratObject::GetAssayData(sce, assay = assay, slot = slot)


  if(slot == "scale.data"|is.null(scale)){
    scale<-FALSE
  }else{
    scale<-scale
  }
  es.max = sctype_score(scRNAseqData = scRNAseqData,
                        scaled = !scale,
                        gs_sctyper = gs_sctyper_list$gs_sctyper_positive,
                        gs_sctyper2 = gs_sctyper_list$gs_sctyper_negative,
                        gene_names_to_uppercase = gene_names_to_uppercase)

  # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
  # In case Seurat is used, it is either sce[["RNA"]]@scale.data (default), sce[["SCT"]]@scale.data, in case sctransform is used for normalization,
  # or sce[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

  # merge by cluster sce@meta.data[, cluster]
  cL_resutls = do.call("rbind", lapply(unique(sce@meta.data[, cluster]), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(sce@meta.data[sce@meta.data[, cluster]==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sce@meta.data[, cluster]==cl)), 10)
  }))


  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
  print(sctype_scores[, 1:3])

  #Please note that sctype_score function (used above) accepts both positive and negative markers through gs_sctyper and gs_sctyper2 arguments.
  #In case, there are no negative markers (i.e. markers providing evidence against a
  #cell being of specific cell type) just set gs_sctyper2 argument to NULL (i.e. gs_sctyper2 = NULL).

  #We can also overlay the identified cell types on UMAP plot:
  #sce@meta.data$sc_typer = ""

  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,]
    sce@meta.data$sc_typer[sce@meta.data[, cluster] == j] = as.character(cl_type$type[1])
  }

  return(sce)
}
