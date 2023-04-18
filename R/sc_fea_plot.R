



#' single cell features plots
#' 
#' Visualize features for specific cell identity classes. First select features with significant differences from the input features.
#' Then use a variety of plot types to visualize above features for specific cell identity classes, including box plots, umap plots, tsne plots,
#' violin plots, and heatmap.
#' @param sce Seurat object after umap and tsne dimensionality reduction
#' @param group Name of a metadata column to group cells by
#' @param assay Assay to pull data from
#' @param slot Data slot to use, options: scale.data, data, counts
#' @param variables Vector of variables to show in plot
#' @param cols Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param path Path of the output saving directory
#' @param fig.type Format of plot saving, such as pdf and png
#' @param dims Vector of dimensionality reduction methods, such as 'umap' and 'tsne'
#' @param index Prefix of file name for saving
#' @param show_variables Number of variables to show in plot, default is 10
#' @param pt.size Point size of 'FeaturePlot', default is 1
#' @param show_box_pvalue If TRUE, paired wised statistical p-value will be shown in boxplot
#' @param show_label If TRUE, labels of cell clusters will be shown in featurePlot 
#' @param show_plot If TRUE, pheatmap of all selected variables will be shown in the Rstudio
#' @param sub_group Name of a metadata column to re-group a subset of cells by
#' @param target Vector of specific cell identity classes in `group` 
#' @param split_by Name of a metadata column to split plot by
#' @param remove_other_celltypes Default is False
#' @param min_cell_count Minimal cell counts in each cluster, default is 10
#' @param dims_for_recluster Number of PCs for re-clustering, default is 30
#' @param palette.heatmap Numeric value corresponding with colors to use for n-colour gradient in heatmap, default is 1, other options: 2, 3, 4, 5, 6
#'
<<<<<<< Updated upstream
=======
#' Visualize features for specific cell identity classes. First select features with significant differences from the input features.
#' Then use a variety of plot types to visualize above features for specific cell identity classes, including box plots, umap plots, tsne plots,
#' violin plots, and heatmap.
#' @param sce Seurat object after umap and tsne dimensionality reduction
#' @param group Name of a metadata column to group cells by
#' @param assay Assay to pull data from
#' @param slot Data slot to use, options: scale.data, data, counts
#' @param variables Vector of variables to show in plot
#' @param cols Vector of colors, users can define the cols manually.  This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param path Path of the output saving directory
#' @param fig.type Format of plot saving, such as pdf and png
#' @param dims Vector of dimensionality reduction methods, such as 'umap' and 'tsne'
#' @param index Prefix of file name for saving
#' @param show_variables Number of variables to show in plot, default is 10
#' @param pt.size Point size of 'FeaturePlot', default is 1
#' @param show_box_pvalue If TRUE, paired wised statistical p-value will be shown in boxplot
#' @param show_label If TRUE, labels of cell clusters will be shown in featurePlot
#' @param show_plot If TRUE, pheatmap of all selected variables will be shown in the Rstudio
#' @param sub_group Name of a metadata column to re-group a subset of cells by
#' @param target Vector of specific cell identity classes in `group`
#' @param split_by Name of a metadata column to split plot by
#' @param remove_other_celltypes Default is False
#' @param min_cell_count Minimal cell counts in each cluster, default is 10
#' @param dims_for_recluster Number of PCs for re-clustering, default is 30
#' @param palette.heatmap Numeric value corresponding with colors to use for n-colour gradient in heatmap, default is 1, other options: 2, 3, 4, 5, 6
#'
>>>>>>> Stashed changes
#' @return Data frame with cells as rows and features as columns
#' @export
#' @author Dongqiang Zeng
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small
#' # Run UMAP map on first 20 PCs
#' pbmc_small <- RunUMAP(object = pbmc_small, dims = 1:10)
#' #Choose features
#' variables<-c("PPBP", "IGLL5", "VDAC3", "CD1C", "AKR1C3", "PF4", "MYL9", "GNLY", "TREML1", "CA2")
#' res<-sc_fea_plot(sce= pbmc_small, variables= variables, group= "RNA_snn_res.1", assay="RNA")
sc_fea_plot<-function(sce,
                      group           = NULL,
                      sub_group       = NULL,
                      target          = NULL,
                      remove_other_celltypes = FALSE,
                      min_cell_count  = 10,
                      show_label      = TRUE,
                      split_by        = NULL,
                      assay           = NULL,
                      slot            = "scale.data",
                      variables       = NULL,
                      show_variables  = 10,
                      dims            = c("umap", "tsne"),
                      pt.size         = 1,
                      cols            = "normal",
                      palette         = 1,
                      show_cols       = FALSE,
                      palette.heatmap = 1,
                      path            = NULL,
                      index           = NULL,
                      show_box_pvalue = T,
                      fig.type        = "pdf",
                      show_plot       = FALSE,
                      dims_for_recluster = 30){
  
  
  if(length(unique(sce@meta.data[, group]))==1) stop("Cell identity has only one level...")
  
  #define-path
  if(is.null(path)){
    file_name<-creat_folder(paste0("0-",group,"-sig-plot"))
  }else{
    file_name<-creat_folder(paste0(path))
  }
  
  #############################################
  cols<- get_cols(cols = cols, palette = palette, show_col = show_cols)
  
  if(is.null(index)){
    prefix<-NULL
  }else{
    prefix<-paste0(index,"-")
  }
  #############################################
  
  if(!is.null(sub_group)) {
    
    
    Idents(sce)<- group
    sces_sub<-subset(sce, idents = target)
    
    if(remove_other_celltypes){
      sces_sub@meta.data[, sub_group]<-ifelse(grepl(sces_sub@meta.data[, sub_group],pattern = target), sces_sub@meta.data[, sub_group], "unassigned" )
    }else{
      sces_sub@meta.data[, sub_group]<-ifelse(grepl(sces_sub@meta.data[, sub_group],pattern = target), sces_sub@meta.data[, sub_group], "other_celltypes" )
    }
    
    #filter out the clusters with few cells
    ###########################################
    sces_sub<-unassign_cell(sce                = sces_sub,
                            cluster            = sub_group,
                            ignore_cell_prefix = NULL,
                            min_cell_count     = min_cell_count,
                            new_col            = NULL,
                            delete_unassigned  = T,
                            return_meta_data   = FALSE)
    ###########################################
    
    sce<-sces_sub
    group<- sub_group
    Idents(sce)<- sub_group
    
    sce <- ScaleData(sce)
    sce <- FindVariableFeatures(object = sce)
    sce <- RunPCA(object = sce, npcs = dims_for_recluster, verbose = TRUE)
    
    set.seed(123)
    sce <- RunTSNE(object = sce, dims = seq(dims_for_recluster), do.fast = TRUE, verbose= T, check_duplicates = FALSE)
    sce <- RunUMAP(sce, reduction = "pca", dims = seq(dims_for_recluster), do.fast = TRUE, verbose= T)
    #########################################
    #########################################
    #
    # p1<-DimPlot(sce, reduction = "tsne", cols = cols, pt.size = pt.size, label = T)
    # p2<-DimPlot(sce, reduction = "umap", cols = cols, pt.size = pt.size, label = T)
    # p<-p1+p2
    #
    # ggsave(p, filename = paste0("0-",prefix,"-subcluster-dimplot-tsne-umap.pdf"), path = file_name$folder_name, width = 12.5, height = 5)
    
  }
  
  
  ##############################################
  reduction_method<-dims
  #############################################
  
  pp<-list(NULL)
  
  for(ii in 1:length(reduction_method)){
    redu<-reduction_method[ii]
    message(paste0(">>>>--- running ", redu))
    
    #' choose cell identity classes-------------------------
    var<- group
    ppi<-dong_dimplot(sce              = sce,
                      reduction        = redu,
                      groups           = var,
                      split.by         = split_by,
                      label            = T,
                      label.size       = 3,
                      pt.size          = 0.5,
                      cols             = cols,
                      seed             = 4321,
                      show_col         = F,
                      width            = 8,
                      height           = 8,
                      w_index          = 7,
                      w_add            = 2,
                      max_category     = 14,
                      show_plot        = F,
                      path             = file_name$folder_name,
                      index            = paste0(prefix,"0-", ii, "-", redu),
                      legend.position  = "right",
                      legend.direction = "vertical",
                      legend.size      = 9,
                      save_plot        = T,
                      fig.type         = fig.type)
    
    pp[[ii]]<-ppi
    
  }
  pp_com<-pp[[1]]+pp[[2]]
  ggsave(pp_com, filename = paste0(prefix,"0-3-combined-umap-tsne-",var,"-",group,".",fig.type),
         width = 18, height = 8, path = file_name$folder_name )
  #####################################################################
  
  DefaultAssay(sce)<-assay
  ##############################################
  input<-extract_sc_data(sce,
                         vars              = variables,
                         assay             = assay,
                         slot              = slot,
                         combine_meta_data = TRUE)
  
  vars<-variables[variables%in%colnames(input)]
  
  
  # boxplot
  if(!slot=="scale.data") input[,vars]<-scale(input[, vars], scale = T, center = T)
  
  
  print(table(input[, group]))
  box_width <- 3 +  length(unique(sce@meta.data[, group])) * 0.5
  
  if(is.null(split_by)){
    feaPlot_width <- 7.5
  }else{
    feaPlot_width <- 7 + 3.5 * length(unique(sce@meta.data[, split_by]))
  }
  
  # statistical analysis
  if(length(unique(input[, group]))<=2){
    
    level<- unique(input[, group])
    level<-level[order(level)]
    res<-batch_wilcoxon(input, target = group, feature = vars)
    res<-rownames_to_column(res, var = "features")
    writexl::write_xlsx(res, paste0(file_name$abspath, prefix, "0-statistical-res-with-",group,".xlsx"))
  }else{
    
    aa<-lapply(input[,vars], function(x) kruskal.test(x~input[, group]))
    res<-data.frame(p.value = sapply(aa, getElement, name = "p.value"),
                    sig_names = vars,
                    statistic = sapply(aa, getElement, name = "statistic"))
    res$p.adj<-p.adjust(res$p.value,method = "BH",n=length(res$p.value))
    res<-res[order(res$p.adj,decreasing = F),]
    res<-rownames_to_column(res, var = "features")
    writexl::write_xlsx(res, paste0(file_name$abspath, prefix, "0-statistical-res-with-",group,".xlsx"))
    #######################################
  }
  
  if(length(vars) > show_variables){
    show_vars<-res$sig_names[1:show_variables]
  }else{
    show_vars<- res$sig_names
  }
  
  message(">>>>--- Features that will be processed:  ")
  print(show_vars)
  
  ######################################################################################
  for (z in 1:length(show_vars)) {
    
    var<-show_vars[z]
    
    message(paste0(">>> Processing feature:  ", var))
    
    width_box<- 4.5+ length(unique(sce@meta.data[, group]))*0.5
    ###################################
    p1<-IOBR:: sig_box(data = input, signature = var, variable = group, cols = cols,
                       
                       palette = "jama", angle_x_text = 60, hjust = 1, show_pvalue = show_box_pvalue, return_stat_res = FALSE)
    ggsave(p1, filename = paste0(prefix, z,"-7-",var,"-boxplot-with-",group, ".",fig.type),
           width = width_box, height = 11, path = file_name$folder_name)
    #####################################
    p1<-IOBR::sig_box(data = input, signature = var, variable = group, cols = cols,
                      palette = "jama", angle_x_text = 60, hjust = 1, show_pvalue = FALSE, return_stat_res = FALSE)
    ggsave(p1, filename = paste0(prefix, z,"-8-",var,"-boxplot-with-",group, ".",fig.type),
           width = width_box, height = 11, path = file_name$folder_name)
    ####################################
    p1<-IOBR::sig_box(data = input, signature = var, variable = group, cols = cols,
                      palette = "jama", angle_x_text = 60, hjust = 1, show_pvalue = FALSE, return_stat_res = TRUE)
    writexl::write_xlsx(p1, paste0(file_name$abspath, prefix, z, "-9-statistical-res-",var,"-",group,".xlsx"))
    
    if(!is.null(assay)) DefaultAssay(sce)<- assay
    print(paste0("Default assay is ", DefaultAssay(sce) ))
    
    print(paste0(">>>-- Colors could be change by parameter: 'cols'"))
    ###############################################################
    
    Idents(sce)<- group
    
    p<- ggplot(input, aes(x= !!sym(group), y = !!sym(var), fill = !!sym(group) )) + geom_violin(trim=FALSE)+ scale_fill_manual(values = cols ) + design_mytheme()
    # help("VlnPlot")
    # VlnPlot(object = sce, features = var, assay = assay,  cols = cols, ncol = 1) & theme(plot.title = element_text(size = 10))
    ggsave(plot = p,filename=paste0(prefix,z,"-", "4-",var,"-VlnPlot_subcluster_markers.", fig.type),
           width = box_width, height = 8,
           path = file_name$folder_name)
    
    # help("FeaturePlot")
    pp<- FeaturePlot(object = sce, features=var, reduction = "umap", split.by = split_by,  label = show_label, pt.size = pt.size, ncol = 1,
                     cols = c("lightgrey", "darkred")) & theme(plot.title = element_text(size = 10))
    ggsave(plot = pp, filename=paste0(prefix,z,"-", "5-",var,"-FeaturePlot-umap.",fig.type),
           width = feaPlot_width, height = 7,
           path = file_name$folder_name)
    
    pp<- FeaturePlot(object = sce, features=var, reduction = "tsne",  split.by = split_by,   label = show_label, pt.size = pt.size, ncol = 1,
                     cols = c("lightgrey", "darkred"))  & theme(plot.title = element_text(size = 10))
    ggsave(plot = pp, filename=paste0(prefix,z,"-","6-",var,"-FeaturePlot-tsne.",fig.type),
           width = feaPlot_width, height = 7,
           path = file_name$folder_name)
    ################################################################
    
  }
  
  
  pp<-      dong_heatmap(sce             = sce,
                         group           = group,
                         feas            = show_vars,
                         assay           = assay,
                         slot            = slot,
                         cols            = cols,
                         palette         = 1,
                         palette.heatmap = palette.heatmap,
                         prefix          = index,
                         disp.max        = 2.5,
                         disp.min        = -2.5,
                         path            = file_name$folder_name,
                         show_plot       = show_plot,
                         show_col        = show_cols,
                         fig.type        = fig.type)
  # Idents(sce)<- group
  # wwidth<- 4 + length(unique(sce@meta.data[,group]))*0.4
  # hheight<- 3 + length(unique(show_vars)) * 0.3
  # ###################################
  # mapal <- colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256)
  # pp<- DoHeatmap(object =  sce,
  #                features     = show_vars,
  #                angle        = 60,
  #                size         = 3.5,
  #                group.colors = cols,
  #                assay        = assay,
  #                slot         = slot)+
  #   scale_fill_gradientn(colours = rev(mapal))
  #
  # if(show_plot) print(pp)
  # ggsave(plot =  pp, filename=paste0(prefix,"0-4-DoHeatmap-",group,".",fig.type),
  #        path = file_name$folder_name,
  #        width = wwidth, height = hheight)
  ###################################
  message(">>>--- Processing VlnPlot ")
  height_VlnPlot<- 4.0+ length(show_vars)/1.5
  width_vln<- 4.5+ length(unique(sce@meta.data[, group]))*0.6
  ###################################
  p<-VlnPlot(object = sce,
             # assay = "RNA",
             # group.by = "scpred_seurat",
             add.noise = FALSE,
             features = show_vars,
             # log = log,
             stack = T,
             flip = T,
             y.max = 4,
             sort = "increasing",
             cols = cols, #palettes(category = "random", palette = 1),
             ncol = 4)& theme(plot.title = element_text(size = 10), legend.position = "none")
  ggsave(p, filename=paste0(index,"-0-6-VlnPlot_subcluster_markers.", fig.type),
         width = width_vln, height = height_VlnPlot,
         path = file_name$folder_name)
  
  ###################################
  message(">>>--- Processing pheatmap_average ")
  
  height_pheatmap<- 4.0+ length(show_vars)/3
  
  res<-pheatmap_average(sce             = sce,
                        assay                = assay,
                        slot                 = slot,
                        marker_res           = NULL,
                        show_features        = vars,
                        top_n                = 20,
                        group                = group,
                        character_limit      = 20,
                        path                 = file_name$folder_name,
                        cols                 = cols,
                        seed                 = 54321,
                        fig.type             = fig.type,
                        file_name_prefix     = paste0(index, "-0-5"),
                        show_col             = FALSE,
                        width                = 12,
                        height               = height_pheatmap )
  ####################################################
  
  return(input)
}


