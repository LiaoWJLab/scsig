
#' Find and show gene expression markers for all identity classes
#'
#'This function can be used to find markers (differentially expressed genes) for each of the identity classes in a dataset
#'based on `FindAllMarkers(){Seurat}`. Then, use `DoHeatmap(){Seurat}` and `pheatmap_average()`to draw a heatmap
#'and `VlnPlot(){Seurat}` to draw a violin plot of single cell feature expression.
#'Besides, `enrich_cluster_degs()` is used to perform enrichment analysis.
<<<<<<< Updated upstream
=======
#'
>>>>>>> Stashed changes
#' @param sce Seurat object
#' @param group A vector of variables to group cells
#' @param verbose Whether to print a progress bar once expression testing begins
#' @param fig.type Format of plot saving, such as pdf and png
#' @param pt.size Size of point
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character, such as normal and random,  to a palette as specified by `IOBR::palettes`. See `palettes` for details
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param show_col Whether to display the palettes
#' @param show_genes Number of genes displayed in the heatmap drawn by `DoHeatmap()`, default is 10
#' @param hwidth width of plot when saving
#' @param show_plot Whether to display the plot
#' @param path Path of the output saving directory
#' @param hheight Height of plot when saving
#' @param assay Assay to use in differential expression testing
#' @param slot Slot to pull data from; note that if test.use is "negbinom", "poisson", or "DESeq2", slot will be set to "counts"
#' @param logfc.threshold Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is 0.15 Increasing logfc.threshold speeds up the function, but can miss weaker signals
#' @param test.use Denotes which test to use. Available options are: "wilcox"(default), "bimod", "roc", "t", "nr", "negbinom", "poisson", "LR", "MAST", "DESeq2".
#' @param only.pos Only return positive markers, default is T
#' @param show_features Number of genes displayed in the violin plot, default is 8
#' @param character_limit Limit length of character, default is 50
#' @param show_genes_pheatmap Number of genes displayed in the heatmap drawn by `pheatmap_average()`, default is 5
#' @param feature_type Feature type, default is gene, otherwise the enrichment analysis will not be done.
#' @param enrich_cutoff_logfc Cutoff for filtering gene signatures which show on average, less than X-fold difference (log-scale) between the two groups of cells for enrichment analysis
#' @param enrich_cutoff_padj Cutoff for filtering genes which adjust P values are < n
#' @param dims_for_recluster Number of PCs to compute and store for re-clustering
#' @param resolution_for_recluster Value of the resolution parameter for re-clustering
#' @param assay_for_recluster Name of assay to use, such as RNA,SCT,integrated
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param feas Genes to test. Default is to use all genes
<<<<<<< Updated upstream
#' @param group_after_recluster
#' @param adjust_assay
#' @param re_scale_tsne_umap
=======
#' @param group_after_recluster if TRUE, finding cluster  will be proceeded in new data
#' @param adjust_assay default is FALSE
#' @param re_scale_tsne_umap if TRUE, dimension reduction will be proceeded in new data
#' @param recluster default is FALSE
>>>>>>> Stashed changes
#'
#' @return
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data("pbmc_small")
#' dong_find_markers(sce = pbmc_small)
dong_find_markers<-function(sce,
                            assay                    = NULL,   #RNA, sct, integrated
                            adjust_assay             = FALSE,
                            slot                     = "scale.data", # count, scale.data
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
                            show_features            = 8,
                            show_plot                = T,
                            path                     = NULL,
                            character_limit          = 50,
                            re_scale_tsne_umap       = FALSE,
                            recluster                = FALSE,
                            assay_for_recluster      = NULL,
                            dims_for_recluster       = 20,
                            group_after_recluster    ="default",
                            resolution_for_recluster = 0.2){

  # -------------------------------------------------------------------------


  # -------------------------------------------------------------------------

  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("FindMarkers-dong")
  }

  path<-creat_folder(file_store)

  # if(!file.exists(file_store)) dir.create(file_store)
  # abspath<-paste(getwd(),"/",file_store,"/",sep ="" )

  message(">>>---Assay used to estimation:")
  if(!is.null(assay)){
    print(paste0(">>>>> ",assay))
  }else{
    print(DefaultAssay(sce))
    DefaultAssay(sce)<- assay
  }


  # find marker genes of clusters------------------------------------------------------
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

  if(tolower(assay)=="sct"&&adjust_assay) sce<- PrepSCTFindMarkers(sce)
  # help("FindAllMarkers")
  ###################################
  sce.markers <- FindAllMarkers(object          = sce,
                                slot            = slot,
                                assay           = assay,
                                features        = feas,
                                only.pos        = only.pos,
                                min.pct         = 0.25,
                                thresh.use      = 0.25,
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
      if(assay_for_recluster!="integrated")  sce <- NormalizeData(object = sce)
      sce <-  ScaleData(sce)
      sce <- FindVariableFeatures(object = sce)

    }

    # sce <- ScaleData(object = sce)
    sce <- RunPCA(object = sce, npcs = dims_for_recluster, verbose = TRUE)

    #findclusters--------------
    sce <- FindNeighbors(sce, dims = seq(dims_for_recluster))
    sce <- FindClusters(sce, resolution = resolution_for_recluster)


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

    if(slot=="scale.data"){
      log = FALSE
    }else{
      log = TRUE
    }

    #############################################################
    if(!is.null(assay)) DefaultAssay(sce)<- assay
    print(paste0("Default assay is ", DefaultAssay(sce) ))

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
            sort = "increasing",
            cols = mycols, #palettes(category = "random", palette = 1),
            ncol = 4)& theme(plot.title = element_text(size = 10), legend.position = "none")
    ggsave(filename=paste0(i+2,"-0-",var,"-VlnPlot_subcluster_markers.", fig.type),
           width = 7,height = 12,
           path = path$folder_name)
    ################################################################

    VlnPlot(object = sce, features = markers_genes,  log = log , cols = mycols, ncol = 4)& theme(plot.title = element_text(size = 10))
    ggsave(filename=paste0(i+2,"-1-",var,"-VlnPlot_subcluster_markers.", fig.type),
           width = 18,height = 13,
           path = path$folder_name)

    FeaturePlot(object = sce, features=markers_genes, reduction = "umap", pt.size = pt.size,
                ncol = 4, cols = c("lightgrey", "darkred")) & theme(plot.title = element_text(size = 10))
    ggsave(filename=paste0(i+2,"-2-",var,"-FeaturePlot_subcluster_markers-umap.",fig.type),
           width = 20,height = 11,
           path = path$folder_name)

    FeaturePlot(object = sce, features=markers_genes, reduction = "tsne", pt.size = pt.size, ncol = 4, cols = c("lightgrey", "darkred"))  & theme(plot.title = element_text(size = 10))
    ggsave(filename=paste0(i+2,"-3-",var,"-FeaturePlot_subcluster_markers-tsne.",fig.type),
           width = 20,height = 11,
           path = path$folder_name)
  }


}
