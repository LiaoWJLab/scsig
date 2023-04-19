




#' Title
#'
#' @param deg
#' @param cluster
#' @param gene
#' @param n
#' @param avg_log2FC
#'
#' @return
#' @export
#'
#' @examples
format_sig_from_df<-function(deg, cluster = "cluster", gene = "gene", avg_log2FC = "avg_log2FC", n = 100){

  # cluster_<- !sym(cluster)
  # avg_log2FC_ <- !sym(avg_log2FC)
  deg  <- as.data.frame(deg)
  deg  <- deg %>% dplyr:: group_by(cluster) %>%  dplyr::top_n(n, avg_log2FC)
  feas <- split(deg, deg[, cluster])
  feas <- lapply(feas, function(x) as.data.frame(x))
  feas <- lapply(feas, function(x) as.character(x[,gene]))
  return(feas)
}




#' draw confusion matrix
#'
#' @param input
#' @param x
#' @param y
#'
#' @return
#' @export
#'
#' @examples
#'
ggplotConfusionMatrix <- function(input, x, y){


  if(class(input)[1]=="Seurat"){
    input<- input@meta.data
  }
  m <-caret:: confusionMatrix(as.factor(input[, x]), as.factor(input[, y]))

  mytitle <- paste("Accuracy", scales:: percent_format()(m$overall[1]),
                   "Kappa", scales:: percent_format()(m$overall[2]))
  p <-
    ggplot(data = as.data.frame(m$table) ,
           aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = log(Freq)), colour = "white") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    geom_text(aes(x = Reference, y = Prediction, label = Freq)) +
    theme(legend.position = "none") + theme_light()+
    ggtitle(mytitle)+
    xlab(x) + ylab(y)
  return(p)
}



#' Title
#'
#' @param sce
#' @param gs
#' @param method
#' @param assay
#' @param slot
#' @param min.feature
#' @param scale
#' @param cluster
#' @param point.size Size of point, default is 1.5
#' @param cols Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param show_col Whether to show color palettes
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param tissue_type
#' @param gs_sctyper
#' @param db_
#' @param db_path
#' @param cell_type
#' @param cell_subset
#' @param study
#' @param gene_names_to_uppercase
#' @param model_scpred
#' @param merge_seurat_cluster
#' @param mini_cluster
#' @param threshold
#' @param source
#' @param path_model_scpred
#' @param show_plot
#' @param save_plot
#' @param path
#' @param width
#' @param height
#' @param deg
#' @param fig.type
#'
#' @return
#' @export
#'
#' @examples
scsig_cell_anno<-function(sce,
                          gs         = NULL,
                          deg        = NULL,
                          method     = c("pca", "sctyper", "model", "censu"),
                          assay      = NULL,
                          slot       = "scale.data",
                          min.feature= 10,
                          scale      = NULL,
                          cluster    = "seurat_clusters",
                          model_scpred= NULL,
                          merge_seurat_cluster = FALSE,
                          mini_cluster         = 50,
                          threshold            = 0.65,
                          source               = "win",
                          path_model_scpred    = NULL,

                          tissue_type= NULL,
                          gs_sctyper = NULL,
                          db_        = "ScTyperDB-merged.xlsx",
                          db_path    = NULL,
                          cell_type  = "base",
                          cell_subset= NULL,
                          study      = NULL,
                          gene_names_to_uppercase = TRUE,
                          cols       = "normal",
                          palette    = 1,
                          show_col   = FALSE,
                          seed       = 123,
                          point.size = 1.5,
                          reduction  = "umap",
                          show_plot  = TRUE,
                          save_plot  = TRUE,
                          width      = 12,
                          height     = 5,
                          path       = NULL,
                          fig.type   = "pdf"){


  message(">>>---Assay used to estimation:")
  if(!is.null(assay)){
    print(DefaultAssay(sce))
  }else{
    print(paste0(">>>>> ",assay))
  }
  #################################
  if(method=="pca"){

    cat(crayon::green(">>>--PCA score of celltype signatures will be estimate and used to celltype annotation...\n"))
    new_cluster<- "pca_cluster"

    if(!is.null(deg)) gs<- format_sig_from_df(deg = deg, cluster = "cluster", gene = "gene")

    sces  <-irGSEA.score(object         = sce,
                         assay          = assay,
                         slot           = slot,
                         update         = FALSE,
                         ncores         = 4,
                         min.cells      = 3,
                         min.feature    = min.feature,
                         custom         = T,
                         geneset        = gs,
                         geneid         = "symbol",
                         method         = "PCAscore")
    #####################################################

    # help("extract_sc_data")
    # pca<- extract_sc_data(sce = sces, assay = "PCAscore",  combine_meta_data = TRUE)

    # help("GetAssay")
    pca<- GetAssay(object = sces, assay = "PCAscore", slot = "scale.data")
    pca <- pca@scale.data
    pca[1:5, 1:5]
    ######################################################
    # cluster<- "RNA_snn_res.0.4"
    # merge by cluster sce@meta.data[, cluster]
    cL_resutls = do.call("rbind", lapply(unique(sces@meta.data[, cluster]), function(cl){
      pca.cl = sort(rowSums(pca[ ,rownames(sces@meta.data[sces@meta.data[, cluster]==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(pca.cl), scores = pca.cl, ncells = sum(sces@meta.data[, cluster]==cl)), 10)
    }))

    sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
    # set low-confident (low ScType score) clusters to "unknown"
    # sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    print(sctype_scores[,1:3])
    #######################################################
    sces@meta.data$pca_cluster = ""
    for(j in unique(sctype_scores$cluster)){
      cl_type = sctype_scores[sctype_scores$cluster==j,]
      sces@meta.data$pca_cluster[sces@meta.data[, cluster] == j] = as.character(cl_type$type[1])
    }
  }


  if(method=="sctyper"){

    cat(crayon::green(">>>-- sctyper will be used to celltype annotation...\n"))
    cat(crayon::green("Reference: Ianevski, A., Giri, A.K. & Aittokallio, T. Nat Commun 13, 1246 (2022). https://doi.org/10.1038/s41467-022-28803-w \n"))

    new_cluster<- "sc_typer"
    sces<- sctyper_anno(sce        = sce,
                        gs         = gs,
                        tissue_type= tissue_type,
                        cell_type  = cell_type,
                        cell_subset= cell_subset,
                        study      = study,
                        assay      = assay,
                        slot       = slot,
                        scale      = scale,
                        cluster    = clsuter,
                        db_        = db_,
                        db_path    = db_path,
                        gene_names_to_uppercase = gene_names_to_uppercase)

  }
  #####################################################
  if(method == "model"){

    new_cluster<- "scpred_celltype"
    cat(crayon::green(">>>-- A scPred model `model_scpred` will be used to celltype annotation...\n"))
    sces<- scpred_cell_anno(sce                  = sce,
                            sce_ref              = model_scpred,
                            name                 = "scpred_celltype",
                            merge_seurat_cluster = merge_seurat_cluster,
                            name_merge_seurat    = cluster,
                            mini_cluster         = mini_cluster,
                            threshold            = threshold,
                            source               = source,
                            path_ref             = path_model_scpred)

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


  p1<-DimPlot(sces, reduction = reduction, label = TRUE, cols = mycols,
              repel = TRUE, pt.size = point.size, group.by = new_cluster)

  p2<-DimPlot(sces, reduction = reduction, label = TRUE, cols = mycols,
              repel = TRUE, pt.size = point.size, group.by = cluster)
  p<-p1+p2

  #################################
  if(show_plot) print(p)

  if(save_plot){
    if(is.null(path)){
      path<- "./"
    }else{
      path<- creat_folder(path)
      path<- path$abspath
    }

    ggsave(p, filename = paste0("1-Celltype-predicted-by-", method, "-", reduction, ".", fig.type), width = width, height = height, path = path)
  }

  # get cell-type by cell matrix
  # scRNAseqData<- SeuratObject::GetAssayData(sce, assay = assay, slot = slot)
  #
  # if(slot == "scale.data"|is.null(scale)){
  #   scale<-FALSE
  # }else{
  #   scale<-scale
  # }
  ##################################


  return(sces)

}
