



#' Modified dimensional reduction plot
#'
#' @param sce Seurat object
#' @param reduction Which dimensionality reduction to use. If NULL, first searches for umap, then tsne, then pca.
#' @param groups Name of one or more metadata columns to group (color) cells by (for example, orig.ident); pass 'ident' to group by identity class
#' @param split.by Name of a metadata column to split plot by
#' @param label Whether to label the clusters
#' @param label.size Size of label
#' @param pt.size Size of point
#' @param cols Vector of colors, each color corresponds to an identity class. This may also be a single character, such as normal and random, to a palette as specified by `IOBR::palettes`. See `palettes` for details
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param show_col Whether to display the palettes
#' @param width Width of plot when saving
#' @param height Height of plot when saving
#' @param w_index Numeric value corresponding to `width`. The parameter works when `split.by` not NULL
#' @param w_add Numeric value corresponding to `width`. The parameter can be used to increase `width` when necessary
#' @param max_category Maximal number of groups
#' @param show_plot Whether to display the plot
#' @param show_plot Whether to display the plot
#' @param path Path of the output saving directory
#' @param index Index number of folder name
#' @param legend.position Position of legend
#' @param legend.direction Direction of legend
#' @param legend.size Font size of legend
#' @param plot_title_size Title size
#' @param axis_title_size Axis title size
#' @param axis_text_size Axis text size
#' @param axis_angle Axis angle
#' @param hjust Horizontal justification
#' @param theme Theme of plot.
#' @param save_plot Whether to save the plot
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param fig.type Format of plot saving, such as pdf and png
#'
#' @return A list of ggplot objects
#' @author Dongqiang Zeng
#' @export
#'
#' @examples
#' data("pbmc_small")
#' scsig_dimplot(object = pbmc_small)
scsig_dimplot<-function(sce,
                       reduction        = "umap",
                       groups           = "orig.ident",
                       split.by         = NULL,
                       label            = T,
                       label.size       = 5,
                       pt.size          = 0.5,
                       cols             = "normal",
                       seed             = 123,
                       palette          = 1,
                       show_col         = FALSE,
                       width            = 9.2,
                       height           = 8.2,
                       w_index          = 7,
                       w_add            = 2,
                       max_category     = 18,
                       show_plot        = T,
                       path             = NULL,
                       index            = 1,
                       theme            = "classic",
                       legend.position  = "right",
                       legend.direction = "horizontal",
                       legend.size      = 0.25,
                       plot_title_size  = 2,
                       axis_title_size  = 2,
                       axis_text_size   = 10,
                       axis_angle       = 0,
                       hjust            = 0.5,
                       save_plot        = TRUE,
                       fig.type         = "pdf"){


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("Dimplot-dong")
  }

  path<-creat_folder(file_store)
  #if(!file.exists(file_store)) dir.create(file_store)
  #abspath<-paste(getwd(),"/",file_store,"/",sep ="" )


  num_cols<-length(unique(as.character(sce@meta.data[,groups])))
  ##############################################
  if(length(cols)==1){
    if(cols=="random"){

      mycols<-palettes(category = "random", palette = palette, show_col = show_col)

      if(length(mycols)<num_cols){
        mycols<-c(mycols, palettes(category = "random", palette = 3, show_col = show_col))
      }

      mycols<-mycols[1:num_cols]
      message(">>>> Default seed is 123, you can change it by `seed`(parameter).")
      set.seed(seed)
      mycols<-mycols[sample(length(mycols), length(mycols))]
      ##########################################

      if(show_col) scales::show_col(mycols)

    }else if(cols == "normal"){

      mycols<-palettes(category = "random", palette = palette, show_col = show_col)
      if(length(mycols)<num_cols){
        mycols<-c(mycols, palettes(category = "random", palette = 3, show_col = show_col))
      }

      mycols<-mycols[1:num_cols]
    }
  }else{
    mycols<-cols
    if(show_col) scales::show_col(mycols)
  }
  ################################################


  for (i in 1:length(groups)) {

    group1<-groups[i]
    message(paste0(">>> Processing group:: ", group1))

    if(!group1%in%colnames(sce@meta.data)) stop(">>>>> groups not exist in meta.data of seurat object!")

    pp <- Seurat:: DimPlot(sce,
                           reduction  = reduction,
                           group.by   = group1,
                           label      = label,
                           split.by   = split.by,
                           label.size = label.size,
                           pt.size    = pt.size,
                           cols       = mycols)

    if(dim(sce)[2]>30000){
      limitsize<-TRUE
    }else{
      limitsize<-FALSE
    }

    if(length(unique(sce@meta.data[,group1]))> max_category & legend.position!="bottom"){

      w_add<-round(length(unique(sce@meta.data[,group1]))/max_category, 0) * w_add
      width<- width + w_add

    }

    if(!is.null(split.by)){
      width<-length(unique(sce@meta.data[,split.by]))*w_index
      fig.name<-paste0("-split-by-",split.by)

      if(length(unique(sce@meta.data[,split.by]))>3){
        limitsize<-FALSE
      }else{
        limitsize<-TRUE
      }
    }else{
      fig.name<-NULL
    }

    if(is.null(index)) index<-1
    #####################################

    mytheme<-design_mytheme(
      theme            = theme,
      plot_title_size  = plot_title_size,
      axis_title_size  = axis_title_size,
      axis_text_size   = axis_text_size,
      axis_angle       = axis_angle,
      hjust            = hjust,
      legend.position  = legend.position,
      legend.direction = legend.direction,
      legend.size      = legend.size
    )
    pp<-pp+mytheme #+guides(color = guide_legend(override.aes = list(size = legend.size)), size = "none")

    if(show_plot) print(pp)
    ####################################


    if(legend.position=="bottom") height<- height + 6

    if(save_plot){
      ggsave(pp,
             filename = paste0(index,"-",i,"-",reduction,"-",group1, fig.name,'.',fig.type),
             path = file_store,
             width = width,
             height = height,
             dpi = 300,
             limitsize = limitsize)

      message(paste0(">>>>>> Figure name is:: ", paste0(index,"-",i,"-",reduction,"-",group1, fig.name,'.',fig.type)))

    }

  }

  return(pp)
}



