


#' Cell type annotation by scType
#'
#' Cell type annotation using specific marker combinations from single-cell transcriptomic data
#' refer to [ScType](https://www.nature.com/articles/s41467-022-28803-w)
#'
#' @param sce Seurat object
#' @param assay Assay to pull from, e.g. RNA, SCT, integrated
#' @param slot Data slot to use, choose from 'counts', 'data', or 'scale.data'
#' @param scale Whether the matrix is scaled，default is NULL
#' @param cluster A vector of variables to group cells by
#' @param db_ Database of manually collected cell type annotation, default is "ScTypeDB_full.xlsx" deposited in data
#' @param tissue_type Default is NULL
#' @param cell_type Cell types options: "base", "epithelial", "myeloid", "tcell", "bcell", "fibroblast" and "endothelial", default is "base"
#' @param cell_subset Default is NULL
#' @param gs_list_positive gene set list that up regulated in cell types
#' @param gs_list_netative gene set list that down regulated in cell types
#' @param gs_data gene signature data with data frame format
#' @param study an option to choose gene signatures
#' @param db_path path of signature data, an example: paste0(base::system.file("data", package = "scsig"),"/ScTyperDB-merged.xlsx")
#' @param gene_names_to_uppercase if TRUE, all gene symbol will be in uppercase
#'
#' @return Seurat object with updated metadata containing cell type annotation
#' @export
#'
#' @examples
#' tnbc<-sctyper_anno(sce = sce, tissue = "Immune system")
sctyper_anno<-function(sce,
                       gs_list_positive        = NULL,
                       gs_list_netative        = NULL,
                       gs_data                 = NULL,
                       tissue_type             = NULL,
                       cell_type               = "base",
                       cell_subset             = NULL,
                       study                   = NULL,
                       assay                   = NULL,
                       slot                    = "scale.data",
                       scale                   = NULL,
                       cluster                 = "seurat_clusters",
                       db_                     = "ScTyperDB-merged.xlsx",
                       db_path                 = NULL,
                       gene_names_to_uppercase = TRUE){


  # prepare gene sets

  if(!is.null(gs_list_positive)|!is.null(gs_list_netative)){

    gsp<- gs_list_positive
    gsn<- gs_list_netative

  }else{

    if(!is.null(db_path)){
      db_<- db_path
    }else{
      db_<- paste0(base::system.file("data", package = "scsig"),"/", db_)
    }

    gs_sctyper_list = gene_sets_prepare(gs_data = gs_data, path_to_db_file = db_, cell_type = cell_type, tissue_type = tissue_type, cell_subset = cell_subset, study = study)

    gsp <- gs_sctyper_list$gs_positive
    gsn <-  gs_sctyper_list$gs_negative

  }

  ###################################

  message(">>>---Assay used to estimation:")
  if(is.null(assay)){
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

  ###################################

  es.max = sctype_score(scRNAseqData = scRNAseqData,
                        scaled = !scale,
                        gs = gsp,
                        gs2 = gsn,
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
