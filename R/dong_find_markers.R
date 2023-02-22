







#' dong_find_markers
#'
#' @param sce
#' @param group
#' @param verbose
#' @param fig.type
#' @param pt.size
#' @param cols
#' @param seed
#' @param show_col
#' @param show_genes
#' @param hwidth
#' @param show_plot
#' @param path
#' @param hheight
#' @param assay
#' @param slot
#' @param logfc.threshold
#' @param test.use
#' @param only.pos
#' @param show_features
#' @param character_limit
#' @param show_genes_pheatmap
#' @param feature_type
#' @param enrich_cutoff_logfc
#' @param enrich_cutoff_padj
#' @param recluster
#' @param dims_for_recluster
#' @param resolution_for_recluster
#' @param assay_for_recluster
#' @param palette
#' @param feas
#' @param group_after_recluster
#' @param adjust_assay
#' @param re_scale_tsne_umap
#'
#' @return
#' @export
#'
#' @examples
dong_find_markers<-function(sce,
                            assay                    = NULL,   #RNA, sct, integrated
                            adjust_assay             = FALSE,
                            slot                     = "scale.data", # count = data, scale.data,
                            feas                     = NULL,
                            group                    = NULL,
                            verbose                  = FALSE,
                            fig.type                 = "pdf",
                            pt.size                  = 1,
                            cols                     = "normal",
                            seed                     = 123,
                            show_col                 = F,
                            palette                  = 1,
                            show_genes               = 10,
                            show_genes_pheatmap      = 5,
                            logfc.threshold          = 0.15,
                            test.use                 = "wilcox",
                            only.pos                 = TRUE,
                            feature_type             = "gene",
                            enrich_cutoff_logfc      = 0.2,
                            enrich_cutoff_padj       = 0.05,
                            hwidth                   = 10,
                            hheight                  = NULL,
                            show_features            = 8,   # print
                            show_plot                = T,
                            path                     = NULL,
                            character_limit          = 50,
                            re_scale_tsne_umap       = FALSE,
                            recluster                = FALSE,
                            assay_for_recluster      = NULL,
                            dims_for_recluster       = 30,
                            group_after_recluster    ="default",
                            resolution_for_recluster = 0.2,
                            slot_vln                 = "scale.data",
                            slot_fea                 = "data"){

  # -------------------------------------------------------------------------


  # -------------------------------------------------------------------------

  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("res-FindMarkers")
  }

  path<-creat_folder(file_store)

  # if(!file.exists(file_store)) dir.create(file_store)
  # abspath<-paste(getwd(),"/",file_store,"/",sep ="" )

  message(">>>---Assay used to estimation:")

  if(!is.null(assay)){
    DefaultAssay(sce)<- assay
    print(paste0(">>>>> ", DefaultAssay(sce)))
  }else{
    print(paste0(">>>>> ", DefaultAssay(sce)))
  }


  # 寻找与cluster之间的marker基因------------------------------------------------------
  if(!is.null(group)){
    message(paste0(">>> Idents of Seurat object is: ", group))
    Idents(sce) <- group
    print(table(as.character(Idents(sce))))
    group2<-group
  }else{
    group2<-"default"
    message(paste0(">>> Idents of Seurat object is: Default... \n"))
    print(table(as.character(Idents(sce))))
  }

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

  #remove features with large name
  ##################################
  # help("PrepSCTFindMarkers")

  if(tolower(assay)=="sct" & adjust_assay) sce<- PrepSCTFindMarkers(sce)
  # help("FindAllMarkers")
  ###################################
  sce.markers <- FindAllMarkers(object          = sce,
                                slot            = slot,
                                assay           = assay,
                                features        = feas,
                                only.pos        = only.pos,
                                min.pct         = 0.15,
                                thresh.use      = 0.15,
                                verbose         = verbose,
                                logfc.threshold = logfc.threshold,
                                test.use        = test.use,
                                fc.name         = "avg_log2FC")

  # DT::datatable(sce.markers)

  print(head(sce.markers))

  prefix<-gsub(group2,pattern = "\\.",replacement = "_")
  writexl::write_xlsx(sce.markers, path = paste0(path$abspath,"1-marker-genes-of-",prefix,".xlsx"))

  if(is.null(show_genes)){

    if( length(unique(Idents(sce))) >= 2)  show_genes<- 10
    if( length(unique(Idents(sce))) >= 3)  show_genes<- 9
    if( length(unique(Idents(sce))) >= 4)  show_genes<- 8
    if( length(unique(Idents(sce))) >= 6)  show_genes<- 6

  }

  ####################################
  top10 <- sce.markers %>% dplyr:: group_by(cluster) %>%  dplyr::top_n(show_genes, avg_log2FC)
  ####################################

  if(is.null(hheight)){
    # if(is.null(group)) stop("group must be define")
    hheight<- 4.8 + length(unique(Idents(sce)))*show_genes/7
  }

  ###################################
  mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
  pp<-DoHeatmap(sce,
                features     = top10$gene,
                angle        = 60,
                size         = 3.5,
                group.colors = mycols,

                assay        = assay,
                slot         = slot)+
    scale_fill_gradientn(colours = rev(mapal))

  if(show_plot) print(pp)
  ggsave(pp,filename=paste0("2-markers-heatmap-of-",prefix,".",fig.type),
         path = path$folder_name,
         width = hwidth,
         height = hheight, dpi = 300)
  ###################################


  if(is.null(show_genes_pheatmap)) {
    show_genes_pheatmap<- show_genes
  }

  height_pheatmap<- 4.0 + length(unique(Idents(sce)))*show_genes*0.3


  pheatmap_average(sce             = sce,
                   assay           = assay,
                   slot            = slot,
                   marker_res      = sce.markers,
                   top_n           = show_genes_pheatmap,
                   group           = group,
                   character_limit = character_limit,
                   path            = path$folder_name,
                   cols            = mycols,
                   seed            = 123,
                   fig.type        = fig.type,
                   show_col        = FALSE,
                   width           = 10,
                   height          = height_pheatmap )
  ####################################################

  #if variables is gene, enrichment analysis will be performed
  if(feature_type=="gene"){

    enrich_cluster_degs(sce.markers       = sce.markers,
                        path              = path$folder_name,
                        index             = 88,
                        cutoff_foldchange = enrich_cutoff_logfc,
                        cutoff_p_adj      = enrich_cutoff_padj)
  }

  #####################################################

  if(recluster){
    #######################################################
    message("----->>>--- reclustering cells will be proceeded....----------------------")

    message(paste0(">>>----Assay that will be used to reclustering: ", assay_for_recluster))

    if(is.null(assay_for_recluster)) stop(">>>--- assay_for_recluster must be specified...")

    DefaultAssay(sce)<-assay_for_recluster

    if(assay_for_recluster=="SCT"){
      "The reclustering result is not convinsible!!! Be cautious....."
      sce<-SCTransform(sce, verbose = TRUE)
    }else{
      if(assay_for_recluster!="integrated") sce <- NormalizeData(object = sce)
      sce <-  ScaleData(sce)
      sce <- FindVariableFeatures(object = sce)

    }

    # sce <- ScaleData(object = sce)
    sce <- RunPCA(object = sce, npcs = dims_for_recluster, verbose = TRUE)
    #findclusters--------------
    sce <- FindNeighbors(sce, dims = seq(dims_for_recluster))
    sce <- FindClusters(sce, resolution = resolution_for_recluster)

    #进一步降维
    ########################################
    if(group_after_recluster=="default") Idents(sce) <- group

    set.seed(123)
    sce <- RunTSNE(object = sce, dims = seq(dims_for_recluster), do.fast = TRUE, verbose= T, check_duplicates = FALSE)
    sce <- RunUMAP(sce, reduction = "pca", dims = seq(dims_for_recluster), do.fast = TRUE, verbose= T)
    #########################################
    #########################################

    p1<-DimPlot(sce, reduction = "tsne", cols = mycols, pt.size = pt.size, label = T)
    p2<-DimPlot(sce, reduction = "umap", cols = mycols, pt.size = pt.size, label = T)
    p<-p1+p2


    if(group_after_recluster=="default"){
      width_dim<- 13
    }else{
      width_dim<- 11
    }
    ggsave(p, filename = paste0("0-",prefix,"-subcluster-dimplot-tsne-umap.",fig.type), path = path$folder_name, width = width_dim, height = 5)

  }else{

    if(re_scale_tsne_umap){
      set.seed(123)
      sce <- RunTSNE(object = sce, dims = seq(30), do.fast = TRUE, verbose= T, check_duplicates = FALSE)
      sce <- RunUMAP(sce, reduction = "pca", dims = seq(30), do.fast = TRUE, verbose= T)
    }
    #########################################
    #########################################

    p1<-DimPlot(sce, reduction = "tsne", cols = mycols, pt.size = pt.size, label = T)
    p2<-DimPlot(sce, reduction = "umap", cols = mycols, pt.size = pt.size, label = T)
    p<-p1+p2

    ggsave(p, filename = paste0("0-",prefix,"-subcluster-dimplot-tsne-umap.",fig.type), path = path$folder_name, width = 12.5, height = 5)
  }
  #############################################################

  #' violin-plot-and-feature-plot-of-each-clusters-----------------------------


  vars<- unique(sort(as.character(Idents(sce))))

  for( i in 1:length(vars) ){

    var<-vars[i]
    message(paste0(">>>--- Processing cluster: ", var))
    markers_df <- FindMarkers(object          = sce,
                              assay           = assay,
                              ident.1         = var,
                              features        = feas,
                              slot            = slot,
                              logfc.threshold = logfc.threshold,
                              min.pct         = 0.1,
                              verbose         = verbose)
    print(x = head(markers_df))
    markers_genes =  rownames(head(x = markers_df, n = show_features))  #show_features一般不能改变，因为涉及到下面的排版问题


    ###########################################
    if(is.null(slot_vln)) slot_vln <-  slot
    if(slot_vln=="scale.data"){
      log = FALSE
    }else{
      slot_vln=="data"
      log = TRUE
    }
    #############################################################
    # if(!is.null(assay)) DefaultAssay(sce)<- assay
    print(paste0("Default assay is ", DefaultAssay(sce)))


    print(paste0(">>>-- Colors could be change by parameter: 'cols'"))
    ###############################################################

    VlnPlot(object = sce,
            # assay = "RNA",
            # group.by = "scpred_seurat",
            add.noise = FALSE,
            features = markers_genes,
            # log = log,
            stack = T,
            flip = T,
            slot = slot_vln,
            sort = "increasing",
            cols = mycols, #palettes(category = "random", palette = 1),
            ncol = 4)& theme(plot.title = element_text(size = 10), legend.position = "none")
    ggsave(filename=paste0(i+2,"-0-",var,"-VlnPlot_subcluster_markers.", fig.type),
           width = 7,height = 12,
           path = path$folder_name)
    ################################################################

    VlnPlot(object = sce, features = markers_genes,  log = log,  slot = slot_vln,  cols = mycols, ncol = 4)& theme(plot.title = element_text(size = 10))
    ggsave(filename=paste0(i+2,"-1-",var,"-VlnPlot_subcluster_markers.", fig.type),
           width = 18,height = 13,
           path = path$folder_name)


    ###############################
    if(is.null(slot_fea)) slot_fea <- "scale.data"

    #################################
    FeaturePlot(object = sce, features = markers_genes, reduction = "umap", pt.size = pt.size,
                slot = slot_fea,
                min.cutoff = 0,
                max.cutoff = 20,
                ncol = 4, cols = c("lightgrey", "darkred")) & theme(plot.title = element_text(size = 10))
    ggsave(filename=paste0(i+2,"-2-",var,"-FeaturePlot_subcluster_markers-umap.",fig.type),
           width = 20,height = 11,
           path = path$folder_name)

    FeaturePlot(object = sce, features = markers_genes, reduction = "tsne", pt.size = pt.size,
                slot = slot_fea,
                min.cutoff = 0,
                max.cutoff = 20,
                ncol = 4, cols = c("lightgrey", "darkred"))  & theme(plot.title = element_text(size = 10))
    ggsave(filename=paste0(i+2,"-3-",var,"-FeaturePlot_subcluster_markers-tsne.",fig.type),
           width = 20,height = 11,
           path = path$folder_name)
  }


}
