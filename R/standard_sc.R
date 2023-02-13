


#' standard_Seurat_data_processing
#'
#' @param eset raw data or normalized data
#' @param file_type if file_type = "10X", data_path must be provided.
#' @param data_path absolute path of 10X file
#' @param project project name,
#' @param index number of folder name,
#' @param plot if true, relevant plot will be drawn
#' @param minFeature minimum features of cell
#' @param minCount minimum count of cell
#' @param findmarkers if true, top 6 marker genes of each cluster will be defined and drawn
#' @param nPCs maximum PC of input
#' @param res resolution of t-SNE
#' @param verbose verbose,
#' @param already_normalized logistic vector, if true, `sce` object will not be normalized
#' @param sce_object Seurat object
#' @param vars_dim default is null
#' @param qc_identity default is null
#' @param save_sce_object default is TRUE
#' @param save_path default is null
#' @param species default is human
#' @param cutoff_percent_mt default is 25
#' @param assay default is RNA
#' @param filter_min_cells default is FALSE
#' @param min_cell_count default is 20
#' @param cols default is 'normal'
#' @param palette default is 1, other options: 2,3,4
#' @param show_col default is TRUE, palettes will be displayed
#' @param seed default is 123
#'
#' @return
#' @export
#'
#' @examples
standard_sc<- function(eset               = NULL,
                       sce_object         = NULL,
                       assay              = "RNA",
                       file_type          = "10X",
                       data_path          = NULL,
                       project            = NULL,
                       qc_identity        = NULL,
                       nPCs               = 30,
                       res                = 1.0,
                       verbose            = FALSE,
                       index              = 1,
                       plot               = FALSE,
                       vars_dim           = "Phase",
                       minFeature         = 2000,
                       minCount           = 1000,
                       species            = "human",
                       cutoff_percent_mt  = 25,
                       filter_min_cells   = FALSE,
                       min_cell_count     = 20,
                       findmarkers        = FALSE,
                       cols               = "normal",
                       palette            = 1,
                       show_col           = T,
                       seed               = 123,
                       already_normalized = FALSE,
                       save_path          = NULL,
                       save_sce_object    = TRUE){

  if(is.null(eset) & file_type =="10X"& is.null(sce_object)){
    sce<-Seurat::CreateSeuratObject(Read10X(paste0(data_path)),project)
  }

  if(!is.null(eset)){
    sce <-Seurat::CreateSeuratObject(eset,project)
  }

  if(!is.null(sce_object)){
    sce<- sce_object
  }

  print(sce)

  if(!is.null(save_path)){
    file_name<-creat_folder(save_path)
  }else{
    file_name<-creat_folder(paste0(index,"-",project))
    abspath<-file_name$abspath
  }
  #############################################


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

  message(" >>>-----  Data will be deposite in ",file_name$abspath)
  ##############################################
  path1<-creat_folder(paste0(file_name$folder_name,"/1-QC-and-Dimplot"))
  #############################################
  print(">>>----- Step-1: Quality control")


  if(species == "human"){
    if(length(unique(grepl('^MT-',rownames(sce))))<=1){
      message(paste0("Mitochondrial genes have been removed"))
      mito<-TRUE
    }else{
      mito<-FALSE
      print("Mitochondrial genes: ")
      print(rownames(sce)[grepl('^MT-',rownames(sce))])
    }

    if(length(unique(grepl('^RP[sl]',rownames(sce))))<=1){
      message(paste0("Ribosomal genes have been removed"))
    }else{
      print("Ribosomal genes")
      print(rownames(sce)[grepl("^RPL|^RPS",rownames(sce))])
    }
  }

  sce <- AddMt(sce , species= species )
  #################################################
  # sce[["percent.mt"]] <- Seurat:: PercentageFeatureSet(sce, pattern = "^mt-")
  # rb.genes <- rownames(sce)[grep("^Rp[sl]",rownames(sce))]
  # C<-Seurat:: GetAssayData(object = sce, slot = "counts")
  # percent.ribo <- Matrix::colSums(C[rb.genes,])/Matrix::colSums(C)*100
  # sce <-Seurat:: AddMetaData(sce, percent.ribo, col.name = "percent.ribo")
  # print(head(sce@meta.data))
  ##################################################



  ##################################################
  if(plot){

    sce_plot<-sce

    if(!is.null(qc_identity)){
      # metadata<- sce@meta.data
      Idents(sce_plot)<- sce_plot[[qc_identity]]
      width_index<- length(unique(as.character(sce_plot[[qc_identity]])))
    }else{
      width_index<-NULL
    }

    width_fea<- 8.5
    width_v<-8.5

    if(!is.null(width_index)){
      width_fea<- 6+ width_index*0.4
      width_v<- 6+ width_index*0.2
    }

    plot1 <-Seurat::FeatureScatter(sce_plot, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- Seurat:: FeatureScatter(sce_plot, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    Seurat:: CombinePlots(plots = list(plot1, plot2))
    ggsave(filename = paste0("1-Mitocondrial-genes-and-feature-count.pdf"),
           width = width_fea, height = 5.8,path = path1$folder_name)

    pp<- VlnPlot(sce_plot, features = c("nFeature_RNA", "nCount_RNA", "percent.rp",  "percent.mt"), ncol = 4)
    ggsave(pp,filename = paste0("2-nFeature-and-nCount-ribo-mt-percent.pdf"),
           width = width_v*2, height = 5.8,path = path1$folder_name)

  }

  # Filtering data
  # Idents(sce)<-sce$orig.ident
  # Idents(sce)<-c("nFeature_RNA","nCount_RNA","percent.mt")
  message(">>>------ For cell subset: Default parameters are : minFeature = 2000, minCount = 1000, cutoff_percent_mt = 25")

  if(mito){
    sce <-subset(sce, subset = nFeature_RNA > minFeature & nCount_RNA > minCount)
  }else{
    sce <- subset(sce, subset = nFeature_RNA > minFeature & nCount_RNA > minCount & percent.mt < cutoff_percent_mt)
  }

  message(">>>------ After filtering cells with low features, low count and high expression of mitochondrial genes ")
  print(sce)
  ######################################################

  print(">>>------ Step-2: Data normalization and dimension reduction")

  path2<-path1
  # path2<-creat_folder(paste0(file_name$folder_name,"/2-Normalization-Dimension-reduction"))

  if(!already_normalized){
    # Normalization
    sce <-Seurat:: NormalizeData(sce, normalization.method =  "LogNormalize", scale.factor = 10000, verbose = verbose)
    # sce = NormalizeData(sce,verbose=verbose)
    #######################################################
  }

  message(">>>------ Cell cycle scoring... ")
  # Read in a list of cell cycle markers, from Tirosh et al, 2015.
  # We can segregate this list into markers of G2/M phase and markers of S phase.
  s.genes <- Seurat::cc.genes$s.genes
  s.genes <- s.genes[s.genes %in% rownames(sce)] # genes in dataset
  g2m.genes <- Seurat::cc.genes$g2m.genes
  g2m.genes <- g2m.genes[g2m.genes %in% rownames(sce)] # genes in dataset

  sce <- CellCycleScoring(object = sce, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)
  print(head(sce@meta.data))
  #########################################################


  sce <- FindVariableFeatures(sce,selection.method = "vst", nfeatures = 3000, verbose=verbose)
  sce <- ScaleData(sce, verbose=verbose) #避免nCounts的影响
  sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce), verbose=verbose)

  ##########################################################
  if(plot){
    # mycols<- IOBR:: palettes(category = "random", show_col = F, show_message = F)

    pp1<-ElbowPlot(sce)
    pp2<-DimHeatmap(sce, dims = 1:12, cells = 100, balanced = TRUE)
    ggsave(pp1, filename = paste0("3-1-ElbowPlot.pdf"),width = 5, height = 5.5, path = path2$folder_name)
    ggsave(pp2, filename = paste0("3-2-DimHeatmapPlot.pdf"),width = 12, height = 12, path = path2$folder_name)
  }
  #######################################################
  message("------------------------------------------------------")

  print(">>>------ Step-3: Find clusters")
  #findclusters--------------
  sce <- FindNeighbors(sce, dims = seq(nPCs), reduction = "pca")


  # sce <- FindClusters(sce, resolution = 0.2, verbose = verbose)
  # message(">>>-- Findclusters when resolution = 0.2...")
  # print(summary(as.factor(sce@meta.data$RNA_snn_res.0.2)))
  # message("------------------------------------------------------")
  # sce <- FindClusters(sce, resolution = 0.8,verbose = verbose)
  # message(">>>-- Findclusters when resolution = 0.8...")
  # print(summary(as.factor(sce@meta.data$RNA_snn_res.0.8)))
  # message("------------------------------------------------------")


###########################################################
  for (xx in seq(0.4, 2.0, by = 0.4)) {

    message(paste0("Resolution is ", xx))
    sce <- FindClusters(sce, resolution = xx, verbose = verbose)

    name_cluster<-paste0(assay, "_snn_res.",xx)
    message(paste0(">>>-- Findclusters when resolution = ", name_cluster," ..."))
    print(summary(as.factor(sce@meta.data[,name_cluster])))

  }
  p<-clustree::clustree(sce, alpha = 0.8) +
    theme(legend.position = "bottom") +
    scale_color_brewer(palette = "Set2") +
    ggraph::scale_edge_color_continuous(low = "grey80", high = "red")

  ggsave(p, filename = paste0(3,"-3-finding-best-resolutions.pdf"),
         width = 16, height = 18, path = path2$folder_name )
  #####################################################


  #####################################################
  sce <- FindClusters(sce, resolution = res, verbose = verbose)
  name_cluster<-paste0(assay, "_snn_res.",res)
  message(paste0(">>>-- Findclusters when resolution = ", name_cluster," ..."))
  print(summary(as.factor(sce@meta.data[,name_cluster])))

  if(filter_min_cells){
    sce<-unassign_cell(sce                = sce,
                       cluster            =  name_cluster,
                       ignore_cell_prefix = "null",
                       min_cell_count     = min_cell_count,
                       new_col            = NULL,
                       delete_unassigned  = T,
                       return_meta_data   = FALSE)
  }



  message(">>>-- These results will be used to draw dimension plot....")
  message("------------------------------------------------------")



  #进一步降维
  ########################################
  set.seed(123)
  sce <- RunTSNE(object = sce, dims = seq(nPCs), do.fast = TRUE, verbose=verbose, check_duplicates = FALSE)
  sce <- RunUMAP(sce, reduction = "pca", dims = seq(nPCs), do.fast = TRUE, verbose=verbose)
  #############################################

  # mycols<- IOBR:: palettes(category = "random",show_col = F, show_message = F)
  ########################################
  tsne_pos=Embeddings(sce,'tsne')
  umap_pos=Embeddings(sce,'umap')
  if(plot){
    pp<-DimPlot(sce,
                reduction  = "tsne",
                label      =T,
                label.size = 5,
                pt.size    = 2,
                # split.by ='orig.ident',
                cols       = mycols)
    ggplot2::ggsave(pp,filename = paste0("4-",project,"-tsne-resolution-",res,".pdf"),
                    width = 9.58, height = 8.67, path = path2$folder_name)
  }
  ##############################################
  ##############################################
  phe=data.frame(cell=rownames(sce@meta.data),
                 cluster =sce@meta.data$seurat_clusters)
  table(phe$cluster)
  ###########################################

  head(tsne_pos)
  dat=cbind(tsne_pos,phe,umap_pos)
  head(dat)
 ############################################

  meta<-sce@meta.data
  save(meta,file=paste0(file_name$abspath,"1-",project,'-meta.data.Rdata'))

  if(save_sce_object) save(sce,file = paste0(file_name$abspath,"2-",project,"-sce-Normalized-cluster-data.RData"))
  ##################################################

  if(plot){
    p <- ggplot(dat, aes(x=tSNE_1,y=tSNE_2,color=cluster))+
      geom_point(size=2)+
      scale_color_manual(values = mycols)
    p <- p+stat_ellipse(data=dat,aes(x=tSNE_1,y=tSNE_2,fill=cluster,color=cluster),
                        geom = "polygon",alpha=0.2,level=0.9,type="t",linetype = 2,show.legend = F)+
      coord_fixed()+
      theme_light()
    print(p)
    theme <- theme(panel.grid =element_blank()) +
      theme(panel.border = element_blank(),panel.background = element_blank()) +
      theme(axis.line = element_line(size=1, colour = "black"))
    p=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
    print(p)
    ggplot2::ggsave(filename = paste0('5-',project,"-pretty_tsne-","resolution-",res,".pdf"),
                    width = 9.58, height = 8.67,path = path2$folder_name)

    ################################################

    p <- ggplot(dat, aes(x=UMAP_1,y=UMAP_2,color=cluster))+
      geom_point(size=2)+
      scale_color_manual(values = mycols)
    p <- p+stat_ellipse(data=dat,aes(x=UMAP_1,y=UMAP_2,fill=cluster,color=cluster),
                        geom = "polygon",alpha=0.2,level=0.9,type="t",linetype = 2,show.legend = F)+
      coord_fixed()+
      theme_light()
    print(p)
    theme <- theme(panel.grid =element_blank()) +
      theme(panel.border = element_blank(),panel.background = element_blank()) +
      theme(axis.line = element_line(size=1, colour = "black"))
    p=p+theme+guides(colour = guide_legend(override.aes = list(size=5)))
    print(p)
    ggplot2::ggsave(filename = paste0('6-',project,"-pretty_umap-","resolution-",res,".pdf"),
                    width = 9.58, height = 8.67,path = path2$folder_name)


  }
  ##################################################


  #
  # if(!is.null(tsne_vars)){
  #
  #   for (i in 1:length(tsne_vars)) {
  #
  #     var<-tsne_vars[i]
  #     metadata<-sce@meta.data
  #     Idents(sce) <- metadata[,var]
  #
  #     pp<-DimPlot(sce,
  #                 reduction = "tsne",
  #                 label=T,
  #                 label.size = 5,
  #                 pt.size = 2,
  #                 # split.by = var,
  #                 cols = mycols)
  #     ggplot2::ggsave(pp,filename = paste0( "3-",i+1,"-T-SNE-",var,"-",project,".pdf"),
  #                     width = 11, height = 8.67, path = path2$folder_name)
  #
  #   }
  #
  # }



  if(!is.null(vars_dim)){

    resolution<- res
    #############################################
    #' 查看降维数据与细胞注释的关系--PCA--TSNE---UMAP------
    reduction_method<-c("pca","tsne","umap")
    #############################################

    for (i in 1:length(vars_dim)) {

      var<-vars_dim[i]

      #' 选择需要分析的变量------labels不多的-------------------------
      # groups<-c("orig.ident","final_celltype","tissue_type","seurat_clusters")
      #' 使用双重循环进行-------------------------
      pp<-list(NULL)
      for(jj in 1:length(reduction_method)){

        reduction<-reduction_method[jj]
        message("------------------------------------------------------")
        message(paste0(">>> Processing method:: ", reduction))

        pj<-dong_dimplot(sce              = sce,
                         reduction        = reduction,
                         groups           = var,
                         split.by         = NULL,
                         label            = T,
                         label.size       = 5,
                         pt.size          = 0.5,
                         cols             = mycols,
                         seed             = seed,
                         show_col         = show_col,
                         palette          = palette,
                         width            = 8,
                         height           = 8,
                         w_index          = 7,
                         w_add            = 2,
                         max_category     = 14,
                         show_plot        = T,
                         path             = path2$folder_name,
                         index            = paste0(jj+3,"-resolution-",resolution),
                         legend.position  = "right",
                         legend.direction = "vertical",
                         legend.size      = 0.25,
                         save_plot        = FALSE) #if false, plot will be transport fo pj
        pp[[jj]]<-pj
      }
      pp_com<-pp[[1]]+pp[[2]]+pp[[3]]+patchwork:: plot_layout(nrow = 1, byrow = TRUE)

      ggsave(pp_com, filename = paste0(i+6,"-combine-pca-tsne-umap-of-",var,".pdf"),
             width = 21, height = 6.5, path = path2$folder_name )

    }


  }



  if(findmarkers){

    path3<-creat_folder(paste0(file_name$folder_name,"/3-Marker-genes-of-clusters"))

    print(">>>------ Step-4: Find marker genes of each cluster")


    hheight<- 4.0 + length(unique(Idents(sce)))*6/4
    dong_find_markers(sce        = sce,
                      group      = NULL,
                      assay      = assay,
                      slot       = "scale.data",
                      verbose    = verbose,
                      fig.type   = "pdf",
                      pt.size    = 0.5,
                      cols       = mycols,
                      seed       = seed,
                      show_col   = show_col,
                      palette    = palette,
                      show_genes = 10,
                      hwidth     = 13,
                      hheight    = hheight,
                      show_plot  = T,
                      path       = path3$folder_name)

  }

  # sce = ScaleData(sce,verbose=verbose)
  # sce = FindVariableFeatures(sce,verbose=verbose)
  # sce = RunPCA(sce,verbose=verbose)
  # sce = RunTSNE(sce,dims=seq(nPCs),verbose=verbose)
  # sce = FindNeighbors(sce,dims=seq(nPCs),verbose=verbose)
  # sce = FindClusters(sce,res=res,verbose=verbose)

  return(sce)
}
