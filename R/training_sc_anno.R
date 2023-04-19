



#' Train a prediction model by scPred
#'
#' Given a seurat object with cell-type labels, train a prediction model for cell type annotation by `scPred`.
#' Refer to [scPred](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1862-5)
#' @param sce A seurat object
#' @param group Column in meta.data containing the cell-type labels of each single cell
#' @param remove_cell A vector of variables in `group` you want to remove from the model
#' @param subset_cell Whether to select the subsets of cells as training model according to certain proportion
#' @param proportion Proportion used to select the subsets of cells. The parameter works when subset_cell ="True".
#' @param dims Number of PCs for clustering
#' @param model Classification model supported via caret package. Default is "mda"
#' @param verbose Whether to print a progress bar once expression testing begins
#' @param fig.type Format of plot saving, such as pdf and png
#' @param pt.size Size of point
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character, such as "normal" and "random",  to a palette as specified by `IOBR::palettes`. See `palettes` for details
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param show_col Whether to display the palettes
#' @param show_plot Whether to display the plot
#' @param path Path of the output saving directory
#' @param cores Number of nodes to be forked
#' @param palette Numeric value corresponding with color palette. Default is 4, other options: 1, 2, 3
#' @param assay Assay to use in model training, e.g. RNA, SCT, integrated
#' @param width  width of plot when saving
#' @param height  height of plot when saving
#' @param recorrect Whether to correct for batch effects
#'
#' @return scPred object
#' @export
#'
#' @examples
training_sc_anno<-function(sce,
                           group,
                           assay        = "RNA",
                           recorrect    = FALSE,
                           remove_cell  = NULL,
                           subset_cell  = FALSE,
                           propotion    =  0.2,
                           dims         = 45,
                           model        = "mda",
                           cores        =  1,
                           verbose      = FALSE,
                           fig.type     = "pdf",
                           pt.size      = 1,
                           cols         = "normal",
                           seed         = 123,
                           show_col     = F,
                           palette      = 4,
                           width        = NULL,
                           height       = NULL,
                           show_plot    = T,
                           path         = NULL ){

  if(class(sce) !="Seurat") stop(">>> Input data must be a Seurat object...")


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("training-sc-annotation-by-scPred")
  }

  path<-creat_folder(file_store)

  ##############################################
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
  ################################################



  # 确定训练目标------------------------------------------------------
  if(!is.null(group)){
    target<- group
    message(paste0(">>> Idents of Seurat object is: ", target))
    Idents(sce) <- sce@meta.data[,target]
    print(table(as.factor(sce@meta.data[,target])))
  }


  md<-as.data.frame(sce@meta.data)
  colnames(md)[which(colnames(md)==target)]<-"target"

  md<-md[!is.na(md$target),]
  md<-md[!md$target=="NA",]


  if(!is.null(remove_cell)){

    message(paste0(">>> Removing cell types: /n"))
    message(remove_cell)
    md<-md[!md$target%in%remove_cell, ]
  }



  if(subset_cell){
    input<-random_strata_cells(input= md, group = "target", propotion = propotion)
    cell_input<-rownames(input)

    message(paste0(">>> Subseting cells in each subset randomly: ",propotion*100, "%... "))
    message(paste0(">>> Final training data: "))

    md<-md[rownames(md)%in%cell_input,]
    print(table(as.factor(md$target)))

  }else{

    message(paste0(">>> Final training data: "))
    cell_input<-rownames(md)
    print(table(as.factor(md$target)))

  }

  # print(cell_input[1:10])
  # subset cells
  sce<-subset(sce, cells = as.character(cell_input))
  DefaultAssay(sce) <- assay



  if(tolower(assay)=="sct"){
    sce <- sce
    if(recorrect == TRUE){

      sce <- SplitObject(sce, split.by = "orig.ident")
      sct_features <- SelectIntegrationFeatures(object.list = sce, nfeatures = 2000)
      sces_merged <- merge(sce[[1]], y = sce[2:length(sce)],  project = "project", merge.data = TRUE)
      VariableFeatures(sces_merged) <- sct_features
      pc.num<-50
      #' 利用变异的feature跑PCA
      sces_merged <- RunPCA(object = sces_merged, assay = "SCT", features = sct_features, npcs = pc.num)
      #' harmony会继续利用PCA进行校正
      library(harmony)
      sces_merged <- RunHarmony(object           = sces_merged,
                                assay.use        = "SCT",
                                reduction        = "pca",
                                dims.use         = 1:50,
                                group.by.vars    = "orig.ident",
                                plot_convergence = TRUE)


      sces_merged <- RunUMAP(object = sces_merged, assay = "SCT", reduction = "harmony", dims = 1:pc.num)
      sces_merged <- FindNeighbors(object = sces_merged, assay = "SCT",
                                   reduction = "harmony",
                                   dims = 1:pc.num,
                                   graph.name = "sct")
      sces_merged <- FindClusters(object = sces_merged, resolution = 1, graph.name = "sct")
      sce<-sces_merged

    }
  }else{
    sce <- sce %>%
      NormalizeData()%>%
      FindVariableFeatures(verbose = verbose) %>%
      ScaleData() %>%
      RunPCA(dims = 1:dims, verbose = verbose) %>%
      RunUMAP(dims = 1:dims)
  }



  pp<- DimPlot(sce, group.by = group, reduction = "umap", label = TRUE, cols = mycols, pt.size = pt.size)

  print(pp)

  if(is.null(width)){
    width<- 9
  }else{
    width<- width
  }

  if(is.null(height)){
    height<- 8
  }else{
    height<- height
  }

  ggplot2:: ggsave(pp, filename = paste0("1-Dimplot-umap-",group, ".",fig.type), width = width, height = height, path = path$folder_name)

  pp2<- DimPlot(sce, group.by = "orig.ident", reduction = "umap", label = TRUE,cols = mycols, pt.size = pt.size)
  ggplot2:: ggsave(pp2, filename = paste0("2-Dimplot-umap-","orig.ident", ".",fig.type), width = width, height = height, path = path$folder_name)


  save(sce, file = paste0(path$abspath,"0-",group,"-reference-seurat-data-for-scPred.RData") )
  ####################################

  #' 训练和初步测试预测效果
  sce <- scPred:: getFeatureSpace(sce, group)
  # Secondly, we train the classifiers for each cell using the trainModel function.
  # By default, scPred will use a support vector machine with a radial kernel.
  #' 并行运算
  ##################################

  if(cores>1){
    library(doParallel)
    cl <- makePSOCKcluster(cores)
    registerDoParallel(cl)
    sce <-scPred:: trainModel(sce, model = model, allowParallel = TRUE)
    stopCluster(cl)
  }else{
    sce <-scPred:: trainModel(sce, model = model, allowParallel = TRUE)
  }

  ##################################
  # sce <- trainModel(sce)
  # Training probabilities for each cell in the reference data can be accessed using the get_probabilities method:
  scPred:: get_probabilities(sce) %>% head()
  perfor<- scPred::get_probabilities(sce)%>% head()

  print(perfor)
  if(length(unique(sce[[group]]))>2){
    # perfor<- rownames_to_column(as.data.frame(perfor), var = "id")
    writexl::write_xlsx(perfor, paste0(path$abspath,"1-Performance-of-predict-reference-data.xlsx"))
    p<-scPred::plot_probabilities(sce)
    if(show_plot) print(p)
    ggsave(p, filename = paste0("2-Probability-of-prediction-",group,".pdf"),width = 13,height = 8, path = path$folder_name)
  }

  ###################################
  ###################################
  #' 将预测模型提取出来-这样就可以丢弃sce对象
  scpred <-scPred:: get_scpred(sce)
  save(scpred, file = paste0(path$abspath,"0-",group,"-cell-annotation-model.RData"))

  return(scpred)

}
