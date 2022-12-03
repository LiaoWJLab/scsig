



#' Title
#'
#' @param sce
#' @param group
#' @param remove_cell
#' @param subset_cell
#' @param propotion
#' @param dims
#' @param model
#' @param verbose
#' @param fig.type
#' @param pt.size
#' @param cols
#' @param seed
#' @param show_col
#' @param show_plot
#' @param path
#' @param cores
#' @param palette
#' @param assay
#' @param width
#' @param height
#' @param recorrect
#'
#' @return
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
