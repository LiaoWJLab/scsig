







#' single cell features plots
#'
#' @param sce Seurat object
#' @param group cell identity
#' @param assay assay that stored variables
#' @param slot data type used to estimation, options: scale.data, data, counts
#' @param variables targets that users want to draw plots
#' @param cols specifying colors of cell clusters
#' @param palette default is NULL, options: jama, jco, aaas, set2, and so on, functions derived from IOBR::palettes()
#' @param path default is "0-",group,"-sig-plot", users can specifying the path
#' @param fig.type default is pdf, other options: png
#' @param dims default is 'umap' and 'tsne'
#' @param index default is null
#' @param show_variables Maximum number of variables displayed
#' @param pt.size default is 1, point size of 'FeaturePlot'
#' @param show_box_pvalue default is true, if true, boxplot will show the paired wised statistical p-value
#' @param show_label default is true: featurePlot will show labels of cell clsuters
#' @param show_plot if true, pheatmap of all selected variables will be shown in the R studio
#' @param sub_group
#' @param target
#' @param split_by
#' @param remove_other_celltypes
#' @param min_cell_count
#' @param dims_for_recluster
#' @param palette.heatmap
#'
#' @return
#' @export
#'
#' @examples
#'
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
  cols<- get_cols(cols = cols, palette = palette, show_col = T)

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

    #去除细胞数量很少的clusters
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

    #' 选择需要分析的变量------labels不多的-------------------------
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
                      save_plot        = T)

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
    feaPlot_width <- 7 + 3.5 *length(unique(sce@meta.data[, split_by]))
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

    writexl::write_xlsx(res, paste0(file_name$abspath, prefix, "0-statistical-res-with-",group,".xlsx"))
    #######################################
  }

  if(length(vars)> show_variables){
    show_vars<-res$sig_names[1:show_variables]
  }else{
    show_vars<- res$sig_names
  }

  message(">>>>--- Features that will be proceeded:  ")
  print(show_vars)


  for (z in 1:length(show_vars)) {

    var<-show_vars[z]

    ######################################################################################
    message(paste0(">>> Processing target:  ", var))

    p1<-IOBR:: sig_box(data = input, signature = var, variable = group, cols = cols,

                       palette = "jama", angle_x_text = 60, hjust = 1, show_pvalue = show_box_pvalue, return_stat_res = FALSE)

    ggsave(p1, filename = paste0(prefix, z,"-7-",var,"-boxplot-with-",group, ".",fig.type),
           width = 10, height = 11, path = file_name$folder_name)

    p1<-IOBR::sig_box(data = input, signature = var, variable = group, cols = cols,

                      palette = "jama", angle_x_text = 60, hjust = 1, show_pvalue = FALSE, return_stat_res = TRUE)

    writexl::write_xlsx(p1, paste0(file_name$abspath, prefix, z, "-8-statistical-res-",var,"-",group,".xlsx"))


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
                         path            = file_name$folder_name,
                         show_plot       = show_plot,
                         show_col        = FALSE,
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
  message(">>>--- Processing pheatmap_average ")

  height_pheatmap<- 4.0+ length(show_vars)/3

  pheatmap_average(sce             = sce,
                   assay           = assay,
                   slot            = slot,
                   marker_res      = NULL,
                   show_features   = vars,
                   top_n           = 20,
                   group           = group,
                   character_limit = 20,
                   path            = file_name$folder_name,
                   cols            = cols,
                   seed            = 54321,
                   fig.type        = fig.type,
                   file_name_prefix= paste0(index, "-0-5"),
                   show_col        = FALSE,
                   width           = 12,
                   height          = height_pheatmap )
  ####################################################

  return(input)
}


