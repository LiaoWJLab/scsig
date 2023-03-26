




#' Feature average expression heatmap
#'
#' Draw a heatmap of feature average expression for all identity classes.
#' @param sce Seurat object
#' @param assay Assay to pull from, such as RNA,SCT,integrated
#' @param slot Data slot to use, choose from 'raw.data', 'data', or 'scale.data'
#' @param marker_res A data frame of markers for all identity classes which output from `FindAllMarkers(){Seurat}`
#' @param top_n Number of top markers to show in the heatmap, default is 6
#' @param group A vector of variables to group cells by
#' @param character_limit Limit length of character, default is 30
#' @param path Path of the output saving directory
#' @param seed Seed of the random number generator, default is 54321. The parameter works when cols ="random"
#' @param fig.type Format of plot saving, such as pdf and png
#' @param show_col Whether to show color palettes
#' @param width Width of plot for saving
#' @param height Height of plot for saving
#' @param cols Vector of colors, users can define the cols manually. This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param show_features A list of features you want to displayed in the heatmap
#' @param file_name_prefix Prefix of file name
#' @param palette Numeric value corresponding with colors to use for the color bar. Default is 1, other options: 2, 3, 4
#' @param palette_for_heatmape Numeric value corresponding with colors to use for n-color gradient. Default is 6, other options: 1, 2, 3, 4, 5
#' @param scale.matrix
#'
#' @return A list containing  markers average expression matrix, markers of each cluster, colors of each cluster, average expression heatmap and markers expression tibble.
#' @export
#'
#' @examples
#' data("pbmc_small")
#' # Find markers for all clusters
#' all.markers <- FindAllMarkers(object = pbmc_small)
#' # Draw a heatmap of markers average expression for all clusters
#' res<-pheatmap_average(sce = pbmc_small, assay ="RNA", slot="data", scale.matrix = all.markers, show_col = F)
pheatmap_average<-function(sce,
                           assay,
                           slot,
                           marker_res       = NULL,
                           top_n            = 6,
                           group            = NULL,
                           show_features    = NULL,
                           character_limit  = 30,
                           path             = NULL,
                           cols             = "random",
                           seed             = 54321,
                           show_col         = T,
                           palette          = 1,
                           palette_for_heatmape = 6,
                           scale.matrix     = NULL,
                           fig.type         = "pdf",
                           width            = 13,
                           height           = NULL,
                           file_name_prefix = "0"){


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("Marker-heatmap-average")
  }

  path<-creat_folder(file_store)

  # if(!file.exists(file_store)) dir.create(file_store)
  # abspath<-paste(getwd(),"/",file_store,"/",sep ="" )

  if(!is.null(assay)) DefaultAssay(sce)<- assay

  # find marker genes of clusters------------------------------------------------------
  if(!is.null(group)){
    message(paste0(">>> Idents of Seurat object is: ", group))

    if(!group%in%colnames(sce@meta.data)) stop("group was not exited in colname of seurat object...")
    Idents(sce) <- group
  }else{
    group2<-"default"
    message(paste0(">>> Idents of Seurat object is: Default... \n"))
  }

  print(table(as.character(Idents(sce))))
  id_sce<- Idents(sce)

  ##############################################
  # substitude
  # mycols<- get_cols(cols = cols, palette = palette, show_col = show_cols)
  ##############################################
  if(length(cols)==1){
    if(cols=="random"){

      mycols<-palettes(category = "random", palette = palette, show_col = show_col, show_message = FALSE)
      message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
      set.seed(seed)
      mycols<-mycols[sample(length(mycols), length(mycols))]
      if(show_col) scales::show_col(mycols)

    }else if(cols == "normal"){
      mycols<-palettes(category = "random", palette = palette, show_col = show_col, show_message = FALSE)
    }
  }else{
    mycols<-cols
    if(show_col) scales::show_col(mycols)
  }
  ################################################
  ##########################################################################

  if(is.null(top_n)){

    if(length(unique(Idents(sce)))>=3) top_n<-10
    if(length(unique(Idents(sce)))>=5) top_n<- 8
    if(length(unique(Idents(sce)))>=6) top_n<- 7
    if(length(unique(Idents(sce)))>=7) top_n<- 6
    if(length(unique(Idents(sce)))>=8) top_n<- 5
  }


  if(!is.null(marker_res)){

    marker_res<-marker_res[marker_res$cluster %in% unique(as.character(Idents(sce))), ]

    message(">>>--- Results of DEs..")
    print(head(marker_res))
    # select degs for heatmap
    degs_top5 <- marker_res %>%
      # filter(cluster!= NA) %>%
      dplyr:: group_by(cluster) %>%
      dplyr::top_n(top_n, -p_val_adj) %>%
      dplyr::top_n(top_n, avg_log2FC) %>%
      dplyr:: arrange(cluster, p_val_adj)

    celltypes<-unique(as.character(degs_top5$cluster))
    show_features<- degs_top5$gene
  }else{
    show_features<- show_features
    celltypes<- unique(id_sce)
  }
  print(paste0(">>>>> Features that will be displayed..."))
  print(show_features)
  ###################################################
  avgData<- SeuratObject::GetAssayData(object =  sce, assay = assay, slot = slot)
  avgData<- avgData[show_features, ]

  eset<- as_tibble(avgData)
  print(paste0(">>>>> Head of feature data..."))
  print(head(eset))

  print(paste0(">>>>> Features not exited in matrix data..."))
  # print(summary(show_features%in%rownames(avgData)))
  print(show_features[!show_features%in%rownames(avgData)])

  # Calculate the average  gene expression value of each cluster
  avgData <- avgData %>%
    apply(1, function(x){
      tapply(x, INDEX = id_sce, FUN = mean, na.rm = T) # ExpMean
    }) %>% t
  ##################################################

  # phData <- MinMax(scale(avgData), -3, 3) # z-score
  # phData[is.na(phData)]<-0

  phData<-avgData
  #remove na and INF
  feas<-feature_manipulation(data = phData, feature = colnames(phData))
  phData<-phData[, feas]

  for (x in 1:length(unique(celltypes))) {
    if(length(unique(rownames(phData)))<length(rownames(phData))){
      rownames(phData)[duplicated(rownames(phData))]<-paste0(rownames(phData)[duplicated(rownames(phData))], "_", x)
    }else{
      break
    }
  }
  ###################################################
  for (dd in 1:150) {
    if(length(unique(substring(rownames(phData), 1, character_limit))) < length(rownames(phData))){
      character_limit<-character_limit+1
    }else{
      break
    }
  }
  rownames(phData)<-substring(rownames(phData), 1, character_limit)
  ####################################################
  cluster_colors<- mycols[1:length(celltypes)]
  names(cluster_colors)<-celltypes
  # cluster_colors <- setNames(brewer.pal(8, "Set1"), levels(as.factor(sces$Model1_merge_subcluster)))
  ####################################################

  mapal <- palettes(category = "heatmap", palette = palette_for_heatmape, counts = 200,show_col = show_col)

  if(is.null(height)){
    # if(is.null(group)) stop("group must be define")
    height<- 5 + length(unique(Idents(sce)))*top_n*0.45
  }
  ####################################################
  if(!is.null(marker_res)){
    annotation_row<- data.frame(cluster = degs_top5$cluster, row.names = rownames(phData))
    cluster_pheatmap<-FALSE
  }else{
    cluster_pheatmap<- TRUE
    if(length(celltypes)!=  length(rownames(phData))){
      celltypes2<-celltypes[1:length(rownames(phData))]
      celltypes2[is.na(celltypes2)]<-  rep(celltypes, 20)[1:sum(is.na(celltypes2))]
      annotation_row<- data.frame(cluster = celltypes2, row.names = rownames(phData))
    }else{
      annotation_row<- data.frame(cluster = celltypes, row.names = rownames(phData))
    }
  }
  ######################################################
  if(is.null(scale.matrix)){
    if(slot=="data"){
      scale_phmap<- "row"
    }else{
      scale_phmap<- "none"
    }
  }else if(scale.matrix){
    scale_phmap<- "row"
  }else{
    scale_phmap<- "none"
  }
  # library(pheatmap)
  # pdf(paste0(path$abspath, "2-markers-heatmap-of-average-", group, ".", fig.type), width = width, height = height)
  p<-pheatmap:: pheatmap(
    phData,
    color             = mapal, #colorRampPalette(c("darkblue", "white", "red3"))(99), #heatmap color
    scale             = scale_phmap,
    cluster_rows      = cluster_pheatmap,
    cluster_cols      = T, #clustered by columns
    cellwidth         = 6,
    cellheight        = 6,
    treeheight_col    = 6,
    clustering_method = "complete",
    show_rownames     = T, #show cluster names
    angle_col         = "45",
    annotation_row    = annotation_row,
    annotation_colors = list(cluster = cluster_colors),
    fontsize          = 6.5,
    silent            =  T) #The color for clusters are sames as previous setting
  # dev.off()
  ggsave(p, filename =  paste0(file_name_prefix, "-pheatmap-of-markers-average-", group, ".", fig.type),
         width = width,
         height = height,
         path = path$folder_name)

  res<-list("p_matrix" = phData, "p_anno" = annotation_row, "p_cols" = cluster_colors,  "plot" = p, "eset" = eset)
  return(res)
}


##########################################################################################
# dir("./4-Integrated-data-signature-analysis-mysignatures/4-Epithelial cells-mysignatures")
# (load("./4-Integrated-data-signature-analysis-mysignatures/4-Epithelial cells-mysignatures/0-DE-signatues-collection-of-Epithelial cells.RData"))
# deg_sig<-read_excel("./4-Integrated-data-signature-analysis-mysignatures/4-Epithelial cells-mysignatures/1-Model1_merge_subcluster/1-PCAscore/1-marker-genes-of-Model1_merge_subcluster.xlsx")
# ###########################################################################################
# pheatmap_average(sce = sces,
#                  assay = "PCAscore", slot = "data",
#                  marker_res = deg_sig,
#                  top_n = 5,
#                  group = "Model1_merge_subcluster",
#                  character_limit = 30,
#                  path = "pheatmap",
#                  cols = "random",
#                  seed = 54321,
#                  max_cols = NULL,
#                  fig.type = "pdf",
#                  show_col = FALSE,
#                  width = 13,
#                  height = 15)
##########################################################################################





