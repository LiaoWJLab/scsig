




#' Title
#'
#' @param sce
#' @param assay
#' @param slot
#' @param marker_res
#' @param top_n
#' @param group
#' @param character_limit
#' @param path
#' @param seed
#' @param fig.type
#' @param show_col
#' @param width
#' @param height
#' @param cols
#' @param show_features
#' @param file_name_prefix
#' @param palette
#'
#' @return
#' @export
#'
#' @examples
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

  # 寻找与cluster之间的marker基因------------------------------------------------------
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

  ###########################################
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
    print(marker_res)
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

  print(paste0(">>>>> Head of feature data..."))
  print(head(avgData))

  print(paste0(">>>>> Features not exited in matrix data..."))
  print(summary(show_features%in%rownames(avgData)))
  print(show_features[!show_features%in%rownames(avgData)])


  # 每个基因在每个cluster里的平均值
  avgData <- avgData %>%
    apply(1, function(x){
      tapply(x,INDEX = id_sce, FUN = mean, na.rm = T) # ExpMean
    }) %>% t
  ##################################################
  phData <- MinMax(scale(avgData), -2.5, 2.5) # z-score

  #remove na and INF
  feas<-feature_manipulation(data = phData, feature = colnames(phData))
  phData<-phData[,feas]

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

  mapal <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

  if(is.null(height)){
    # if(is.null(group)) stop("group must be define")
    height<- 5 + length(unique(Idents(sce)))*top_n*0.45
  }

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
  # library(pheatmap)
  # pdf(paste0(path$abspath, "2-markers-heatmap-of-average-", group, ".", fig.type), width = width, height = height)
  p<-pheatmap:: pheatmap(
    phData,
    color             = mapal, #colorRampPalette(c("darkblue", "white", "red3"))(99), #配色
    scale             = "row",
    cluster_rows      = cluster_pheatmap, #不按行聚类
    cluster_cols      = T, #按列聚类
    cellwidth         = 15,
    cellheight        = 15,
    treeheight_col    = 15,
    clustering_method = "complete",
    show_rownames     = T, #显示cluster名
    angle_col         = "45",
    annotation_row    = annotation_row,
    annotation_colors = list(cluster = cluster_colors),
    silent            =  T) #注释的配色就用前面设置的cluster的颜色
  # dev.off()
  ggsave(p, filename =  paste0(file_name_prefix, "-pheatmap-of-markers-average-", group, ".", fig.type),
         width = width, height = height, path = path$folder_name)
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





