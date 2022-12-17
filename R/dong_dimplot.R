



#' Modified_Dimplot
#'
#' @param sce seurat object
#' @param reduction default is umap, other options: tsne, pca
#' @param groups colors of group
#' @param split.by variables used to split the figures
#' @param label true or false
#' @param label.size size of label
#' @param pt.size size of point
#' @param cols users can define the cols manualy
#' @param seed default is 54321
#' @param show_col default is FALSE
#' @param width width of figure
#' @param height height of figure
#' @param w_index
#' @param w_add
#' @param max_category
#' @param show_plot
#' @param path
#' @param index
#' @param legend.position
#' @param legend.direction
#' @param legend.size
#' @param plot_title_size
#' @param axis_title_size
#' @param axis_text_size
#' @param axis_angle
#' @param hjust
#' @param theme
#' @param save_plot
#' @param palette
#'
#' @return
#' @export
#'
#' @examples
dong_dimplot<-function(sce,
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
                       width            = 8,
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
                       save_plot        = TRUE){


  if(!is.null(path)){
    file_store<-path
  }else{
    file_store<-paste0("Dimplot-dong")
  }

  if(!file.exists(file_store)) dir.create(file_store)
  abspath<-paste(getwd(),"/",file_store,"/",sep ="" )


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

    limitsize<-TRUE
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
             filename = paste0(index,"-",i,"-",reduction,"-",group1, fig.name,'.pdf'),
             path = file_store,
             width = width,
             height = height,
             dpi = 300,
             limitsize = limitsize)

      message(paste0(">>>>>> Figure name is:: ", paste0(index,"-",i,"-",reduction,"-",group1, fig.name,'.pdf')))

    }

  }

  return(pp)
}



