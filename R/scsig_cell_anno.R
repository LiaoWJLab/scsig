




#' Gene lists from matrix of gene expression markers
#' Get gene lists from matrix of gene expression markers for all identity classes 
#' @param deg Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
#' @param cluster Name of the column in which the clusters are located
#' @param gene Name of the column in which the markers are located
#' @param n Number of selected top ranked markers
#' @param avg_log2FC Name of the column in which the average log2FC values are located
#'
#' @return A list containing top n gene markers of each cell types
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




#' Create and show a confusion matrix
#' Calculates a cross-tabulation of observed and predicted classes with associated statistics by `confusionMatrix{caret}`,and show it by`ggplot`with 
#' "Accuracy" and "Kappa" values.
#' @param input Seurat object or dataframe 
#' @param x Name of a metadata column as a factor of predicted classes
#' @param y Name of a metadata column as a factor of classes to be used as the true results
#' @param axis_angle Axis angle 
#'
#' @return A plot object 
#' @export
#'
#' @examples
#'
ggplotConfusionMatrix <- function(input, x, y, axis_angle = 60){


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
    # theme(legend.position = "none") + theme_light()+
    design_mytheme(legend.position = "none", axis_angle = axis_angle)+
    ggtitle(mytitle)+
    xlab(x) + ylab(y)
  return(p)
}



#' Single cell annotation by multiple methods
#' 
#' Here, we provide three methods for single cell annotation : [ScType](https://www.nature.com/articles/s41467-022-28803-w) and PCA based on gene markers, 
#' and [scPred](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1862-5) based on reference datasets.
#' @param sce Seurat object
#' @param subcluster Name of main cell type to perform cell subtype annotation, such as "Myeloid, T-cells, B-cells". Default is NULL
#' @param main_celltype Name of the metadata column in which the main cell type annotation is located 
#' @param gs Input own genesets as a list. Each element in the list is a gene set. The parameter works when `method`= "pca" or "sctype"
#' @param method Name of methods, choose from "pca", "sctyper" and  "model"
#' @param assay Assay to pull from, e.g. RNA, SCT, integrated
#' @param slot Data slot to use, choose from 'counts', 'data', or 'scale.data'
#' @param min.feature The minimum genes per cell, default 10. The parameter works if a scRNA-seq matrix is input.
#' @param scale Whether the matrix is scaled，default is NULL
#' @param cluster A vector of variables to group cells by
#' @param point.size Size of point, default is 1.5
#' @param cols Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param show_col Whether to show color palettes
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param reduction Which dimensionality reduction to use. If not specified, first searches for umap, then tsne, then pca
#' @param tissue_type Tissue type. Default is NULL
#' @param gs_data gene signature data with data frame format. The parameter works when method ="sctype" 
#' @param db_ Database of manually collected cell type annotation, default is "ScTypeDB_full.xlsx" deposited in data. The parameter works when method ="sctype"
#' @param db_path Path of gene signatures data, an example: paste0(base::system.file("data", package = "scsig"),"/ScTyperDB-merged.xlsx"). The parameter works when method ="sctype"
#' @param cell_type Cell types options: "base", "epithelial", "myeloid", "tcell", "bcell", "fibroblast" and "endothelial", default is "base". The parameter works when method ="sctype"
#' @param cell_subset Cell subtypes options. Default is NULL. The parameter works when method ="sctype"
#' @param study Studes options to choose gene signatures. The parameter works when method ="sctype"
#' @param gene_names_to_uppercase If TRUE, all gene symbol will be in uppercase.The parameter works when method ="sctype"
#' @param model_scpred A Seurat object with trained model(s) using scPred or a scPred object. The parameter works when method ="model"
#' @param merge_seurat_cluster Whether to merge prediction results by`scPred` with clustering results. The parameter works when method ="model"
#' @param mini_cluster Minimal cell counts of the cluster. The parameter works when method ="model"
#' @param threshold Threshold used for probabilities to classify cells into classes. All cells below this threshold value will be labels as "unassigned". 
#' In the case of binary classification (two cell tyoes), a threshold of 0.5 will force all cells to be classified to any of the two cell types. 
#' For multi-class classification, if there's no probability higher than the threshold associated to a cell type, this will be labelled as "unassigned". The parameter works when method ="model"
#' @param source Character string related to computer system name choose from win and linux. The parameter works when method ="model"
#' @param path_model_scpred The path or the name of the file where the scPred model is read from. The parameter works when method ="model"
#' @param show_plot Whether to show plots
#' @param save_plot whether to save plots
#' @param path Path of the output saving directory
#' @param width Width of plot when saving
#' @param height Height of plot when saving
#' @param deg Matrix containing a ranked list of putative markers, and associated statistics (p-values, ROC score, etc.)
#' @param fig.type Format of plot saving, such as pdf and png
#'
#' @return Seurat object with updated metadata containing cell type annotation
#' @export
#'
#' @examples
scsig_cell_anno<-function(sce,
                          subcluster = NULL,
                          main_celltype = NULL,
                          gs         = NULL,
                          gs_data    = NULL,  # for sctyper
                          deg        = NULL,
                          top_n_deg  = 100,
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
  sce1<-sce
  ###########################
  if(!is.null(subcluster)){
    if ((subcluster %in% sce@meta.data[,main_celltype]) == F) {
      stop(paste0(subcluster, " not in the main celltype."))
    }
    Idents(sce)<-main_celltype
    sce<-subset(sce, idents= subcluster)
  }
  #################################
  if(method=="pca"){

    cat(crayon::green(">>>--PCA score of celltype signatures will be estimate and used to celltype annotation...\n"))
    new_cluster<- "pca_cluster"

    if(!is.null(deg)) gs<- format_sig_from_df(deg = deg, n = top_n_deg)

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

    if(!is.null(deg)) gs<- format_sig_from_df(deg = deg, n = top_n_deg)

    sces<- sctyper_anno(sce                     = sce,
                        gs_list_positive        = gs,
                        gs_data                 = gs_data,
                        tissue_type             = tissue_type,
                        cell_type               = cell_type,
                        cell_subset             = cell_subset,
                        study                   = study,
                        assay                   = assay,
                        slot                    = slot,
                        scale                   = scale,
                        cluster                 = cluster,
                        db_                     = db_,
                        db_path                 = db_path,
                        gene_names_to_uppercase = gene_names_to_uppercase)

  }
  #####################################################
  if(method == "model"){

    new_cluster<- "scpred_celltype_no_rejection"
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
  if(!is.null(subcluster)){
    sce1$sub_celltype<-sce1@meta.data[,main_celltype]
    sce1$sub_celltype[sce1@meta.data[,main_celltype]==subcluster]<-sces@meta.data[,new_cluster]
    sces<-sce1
  }
  #################################
  return(sces)

}
