


#' single cell annotation by scType
#'
#' @param sce seurat object
#' @param assay default is RNA
#' @param slot default is scaled data
#' @param scale default is NULL
#' @param cluster cluster used to annotation
#' @param point.size default is 1.5
#' @param db_ default is ScTypeDB_full.xlsx, deposited in data
#' @param cols 'random', 'normal', or colors
#' @param palette options: 1, 2
#' @param show_col default is TRUE
#' @param seed default is 123 if cols = 'random'
#' @param reduction default is umap
#' @param gs default is null, user can provide gene set file manually
#' @param tissue_type default is null = base
#' @param cell_type default is null,
#' @param cell_subset default is null
#'
#' @return
#' @export
#'
#' @examples
#' tnbc<-sc_type_anno(sce = sce, tissue = "Immune system")
sc_type_anno<-function(sce,
                       gs         = NULL,
                       tissue_type= NULL,
                       cell_type  = "base",
                       cell_subset= NULL,
                       study      = NULL,
                       assay      = NULL,
                       slot       = "scale.data",
                       scale      = NULL,
                       cluster    = "seurat_clusters",
                       cols       = "normal",
                       palette    = 1,
                       show_col   = T,
                       seed       = 123,
                       point.size = 1.5,
                       reduction  = "umap",
                       db_        = "ScTyperDB-merged.xlsx",
                       gene_names_to_uppercase = TRUE){



  db_<- paste0(base::system.file("data", package = "mysc"),"/", db_)

  # prepare gene sets
  gs_list = gene_sets_prepare(gs = gs, path_to_db_file = db_, cell_type = cell_type, tissue_type = tissue_type, cell_subset = cell_subset)
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
                        gs = gs_list$gs_positive,
                        gs2 = gs_list$gs_negative,
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
  print(sctype_scores[,1:3])

  #Please note that sctype_score function (used above) accepts both positive and negative markers through gs and gs2 arguments.
  #In case, there are no negative markers (i.e. markers providing evidence against a
  #cell being of specific cell type) just set gs2 argument to NULL (i.e. gs2 = NULL).

  #We can also overlay the identified cell types on UMAP plot:
  sce@meta.data$sc_typer = ""

  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    sce@meta.data$sc_typer[sce@meta.data[, cluster] == j] = as.character(cl_type$type[1])
  }

  #####################################################
  if(length(cols)==1){
    if(cols=="random"){

      mycols<-palettes(category = "random", palette = palette, show_col = show_col)
      message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
      set.seed(seed)
      mycols<-mycols[sample(length(mycols), length(mycols))]
      if(show_col) scales::show_col(mycols)

    }else if(cols == "normal"){

      mycols<-palettes(category = "random", palette = palette, show_col = show_col)
    }
  }else{
    mycols<-cols
    if(show_col) scales::show_col(mycols)
  }


  p1<-DimPlot(sce, reduction = reduction, label = TRUE, cols = mycols,
              repel = TRUE, pt.size = point.size, group.by = 'sc_typer')

  p2<-DimPlot(sce, reduction = reduction, label = TRUE, cols = mycols,
              repel = TRUE, pt.size = point.size, group.by = cluster)
  p<-p1+p2

  print(p)
  return(sce)
}
