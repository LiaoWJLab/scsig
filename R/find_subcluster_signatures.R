





#' Title
#'
#' @param sce seurat object
#' @param assay default is NULL
#' @param col_celltype name of celltype
#' @param col_sub_celltype name of subcluster
#' @param groups column of
#' @param signature_for_deg signature list
#' @param min_cell_count minimal cell count in each cluster
#' @param no_limitation_celltypes celltypes with no cell count limitation
#' @param sig_methods method used to etimate signature score
#' @param path path to deposided results
#' @param sig_file_name prefix of path
#' @param index default is null, user can choose specific cell types to save computing source
#' @param cols default is normal
#' @param show_feas_pheatmap
#' @param character_limit_heatmap
#' @param remove_other_celltypes
#' @param find_cluster_sig default is FALSE
#' @param groups_for_cluster
#' @param signature_character_limit remove signature with large names
#' @param seed
#' @param palette
#' @param show_col
#'
#' @return
#' @export
#'
#' @examples
find_subcluster_signatures<-function(sce,
                                     assay                     = "integrated",
                                     col_celltype              = "Model1_merge_no_rejection",
                                     col_sub_celltype          = NULL,
                                     groups                    = c("Model1_merge_subcluster", "integrated_snn_res.1"),
                                     cols                      = "normal",
                                     seed                      = 123,
                                     palette                   = 1,
                                     show_col                  = F,
                                     signature_for_deg         = signature_collection,
                                     signature_character_limit = 60,
                                     min_cell_count            = 50,
                                     no_limitation_celltypes   = "Epithelial cells",
                                     index                     = NULL,
                                     sig_methods               = c("PCAscore","ssgsea", "AUCell"),
                                     path                      = "4-Diff-sig-of-subclusters",
                                     sig_file_name             = "mysig",
                                     show_feas_pheatmap        = 8,
                                     character_limit_heatmap   = 50,
                                     remove_other_celltypes    = TRUE,
                                     find_cluster_sig          = FALSE,
                                     groups_for_cluster        = NULL,
                                     recluster                 = FALSE){


  # assay<-"SCT"
  # col_celltype<-"Model1_merge_no_rejection"
  # col_sub_celltype<-"Model1_merge_subcluster"
  # #需要比较的分组
  # groups<-c("Model1_merge_subcluster", "SCT_snn_res.1")
  # signature_for_deg<-signature_collection
  # sig_file_name<- "mysig"
  # sig_methods<-c("PCAscore","ssgsea", "AUCell")

  includ_sig<-nchar(names(signature_for_deg)) <= signature_character_limit
  signature_for_deg<- signature_for_deg[includ_sig]

  #首先构建储藏结果的文件夹和路径
  file_name<-creat_folder(paste0(path, sig_file_name))

  # 4.计算富集分数
  # 当你的ncore设置大于1的时候，发生下面的错误：Error (Valid ‘mctype’: ‘snow’ or ‘doMC’)，你应该检查一下你的AUCell 版本，确保版本大于等于1.14 。如果你比较懒，那你直接把ncore设置为1也是可以的，只是运行速度会稍微慢一点。


  message(">>>---Celltype of Seurat object:")
  print(table(sce[[col_celltype]]))

  celltypes<-unique(sce@meta.data[,col_celltype])


  message(">>>---The order of celltypes:")
  print(celltypes)


  message(">>>--------------------------")

  if(is.null(index)) index<- 1:length(celltypes)

  message(">>>---Celltypes that will be processed:")
  print(celltypes[index])


  message(">>>---Assay used to calculate signature score:")
  if(!is.null(assay)){
    print(DefaultAssay(sce))
  }else{
    print(assay)
  }

  ##################################################################################


  if(find_cluster_sig){

    sce_a<-unassign_cell(sce               = sce,
                        cluster            = col_celltype,
                        ignore_cell_prefix = NULL,
                        min_cell_count     = min_cell_count,
                        new_col            = NULL,
                        delete_unassigned  = T,
                        return_meta_data   = FALSE)

    #############################################
    path_res<-creat_folder(file_name$folder_name, paste0(0,"-","all-clusters","-",sig_file_name))
    #############################################
    ##############################################

    # help("irGSEA.score")
    sce_a <-irGSEA.score(object         = sce_a,
                         assay          = assay,
                         slot           = "scale.data",
                         update         = FALSE,
                         seeds          = 123,
                         ncores         = 4,
                         min.cells      = 3,
                         min.feature    = 0,
                         custom         = T,
                         geneset        = signature_for_deg,
                         msigdb         = T,
                         species        = "Homo sapiens",
                         category       = "H",
                         subcategory    = NULL,
                         geneid         = "symbol",
                         method         = sig_methods,
                         aucell.MaxRank = 2000,
                         ucell.MaxRank  = 2000,
                         kcdf           = 'Gaussian')
    ##########################################

    # 返回一个Seurat对象，富集分数矩阵存放在RNA外的assay中
    Seurat::Assays(sce_a)
    #> [1] "RNA"       "AUCell"    "UCell"     "singscore" "ssgsea"
    # sce_a@assays$ssgsea
    ##########################################

    # 5.整合差异基因集
    # Wlicox test is perform to all enrichment score matrixes and gene sets
    # with adjusted p value &lt; 0.05 are used to integrated through RRA.
    # Among them, Gene sets with p value &lt; 0.05 are statistically
    # significant and common differential in all gene sets enrichment analysis
    # methods. All results are saved in a list.
    result.dge<-list(NULL)
    # groups<-c("Model1_merge_subcluster", "integrated_snn_res.1")


    if(is.null(groups_for_cluster)) groups_for_cluster<- col_celltype

    groups_for_cluster<-groups_for_cluster[groups_for_cluster%in%colnames(sce_a@meta.data)]

    names_group<-NULL
    for (j in 1:length(groups_for_cluster)) {

      group<-groups_for_cluster[j]
      #去掉细胞数量很少的cluster==大于50个
      #########################################
      input<-unassign_cell(sce                = sce_a,
                           cluster            = group,
                           ignore_cell_prefix = NULL,
                           min_cell_count     = min_cell_count,
                           new_col            = NULL,
                           delete_unassigned  = T,
                           return_meta_data   = T)

      #' 如果目标变量没有亚组直接跳过此次循环====================

      print(table(input[, group]))

      if(length(unique(input[, group]))==1) next

      names_group<-c(names_group, group)
      result.dge[[j]] <- irGSEA.integrate(object   = sce_a,
                                          group.by = group,
                                          metadata = NULL,
                                          col.name = NULL,
                                          method   = sig_methods)
    }
    names(result.dge)<-names_group
    class(result.dge)

    # 储存数据
    save(sce_a, result.dge, file = paste0(path_res$abspath, "0-DE-signatues-of-","all-clusters",".RData"))
    ##########################################


    for (jj in 1:length(sig_methods)) {

      DefaultAssay(sce_a) <- sig_methods[jj]

      for(ii in 1:length(groups_for_cluster)){
        group<-groups_for_cluster[ii]

        #去掉细胞数量很少的cluster==大于50个
        #########################################
        sces2<-unassign_cell(sce                = sce_a,
                             cluster            = group,
                             ignore_cell_prefix = NULL,
                             min_cell_count     = min_cell_count,
                             new_col            = NULL,
                             delete_unassigned  = T,
                             return_meta_data   = FALSE)

        #' 如果目标变量没有亚组直接跳过此次循环====================

        print(table(sces2@meta.data[, group]))
        if(length(unique(sces2@meta.data[, group]))==1) next


        path2<-creat_folder(path_res$folder_name, paste0(ii,"-",group), paste0(jj,"-",sig_methods[jj]))



        #############################################
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
        # cols1<-c(palettes(palette="jama"), palettes(palette="set2")[-1])
        #############################################
        p1<-DimPlot(sces2, reduction = "tsne", group.by = group, cols = mycols, pt.size = 0.5, label = T, label.size = 4)
        p2<-DimPlot(sces2, reduction = "umap", group.by = group, cols = mycols, pt.size = 0.5, label = T, label.size = 4)
        p<-p1+p2
        ggsave(p, filename = paste0("0-",group,"-subcluster-dimplot-tsne-umap.pdf"), path = path2$folder_name, width = 13, height = 6 )

        # DefaultAssay(sces) <- sig_methods[jj]
        dong_find_markers(sce                      = sces2,
                          assay                    = sig_methods[jj],
                          slot                     = "scale.data",
                          group                    = group,
                          verbose                  = T,
                          feature_type             = "signature",
                          fig.type                 = "pdf",
                          pt.size                  = 0.5,
                          cols                     = cols,
                          seed                     = seed,
                          palette                  = palette,
                          show_col                 = show_col,
                          show_genes               = 7,
                          show_genes_pheatmap      = show_feas_pheatmap,
                          hwidth                   = 19,
                          hheight                  = NULL,
                          show_plot                = T,
                          path                     = path2$folder_name,
                          character_limit          = character_limit_heatmap,
                          recluster                = recluster,
                          dims_for_recluster       = 8,
                          resolution_for_recluster = 0.3,
                          assay_for_recluster      =  assay)

      }

    }

  }

  #################################################################################
  ###############################++++find_subcluster_signatures####################
  #################################################################################
  #################################################################################
  if(!is.null(col_sub_celltype)){

    for(i in index){

      celltype<-celltypes[i]
      message(paste0(">>>>----Processing celltypes ", celltype))
      ############################################

      #选择细胞大类
      # i<-3
      # celltype<-"B lymphocytes"
      if(!is.null(assay)) DefaultAssay(sce)<-assay

      # print(rownames(sce))

      Idents(sce)<- col_celltype

      sces<-subset(sce, idents = celltype)

      if(remove_other_celltypes){
        sces@meta.data[, col_sub_celltype]<-ifelse(grepl(sces@meta.data[, col_sub_celltype],pattern = celltype), sces@meta.data[, col_sub_celltype], "unassigned" )
      }

      # sces<-subset(sce, Model1_merge_no_rejection == celltype)
      ###########################################

      #去除细胞数量很少的clusters
      ###########################################
      if(!is.null(no_limitation_celltypes)){
        if(celltype == no_limitation_celltypes){
          ignore<- celltype
        }else{
          ignore<- NULL
        }
      }else{
        ignore<- NULL
      }
      ##########################################

      sces<-unassign_cell(sce                = sces,
                          cluster            = col_sub_celltype,
                          ignore_cell_prefix = ignore,
                          min_cell_count     = min_cell_count,
                          new_col            = NULL,
                          delete_unassigned  = T,
                          return_meta_data   = FALSE)

      #' 如果目标变量没有亚组直接跳过此次循环====================
      if(length(unique(sces@meta.data[, col_sub_celltype]))==1) {
        message(">>>--- cells with only one level, this celltype will be skiped...")
        print(table(sces@meta.data[, col_sub_celltype]))
        next
      }

      ############################################
      DefaultAssay(sces)
      #############################################
      path_res<-creat_folder(file_name$folder_name, paste0(i,"-",celltype,"-",sig_file_name))
      #############################################
      # 4.计算富集分数
      # 当你的ncore设置大于1的时候，发生下面的错误：Error (Valid ‘mctype’: ‘snow’ or ‘doMC’)，
      # 你应该检查一下你的AUCell 版本，确保版本大于等于1.14 。如果你比较懒，那你直接把ncore设置为1也是可以的，只是运行速度会稍微慢一点。
      # data("signature_collection")
      # names(signature_collection)[which(names(signature_collection)=="TGF\xa6\xc2.myCAF")]<-"TGFb-myCAF"
      # names(signature_collection)[which(names(signature_collection)=="IFN\xa6\xc3.iCAF")]<-"IFNG-iCAF"
      # signature_collection<-signature_collection[-length(signature_collection)]
      # ##############################################

      # help("irGSEA.score")
      sces <- irGSEA.score(object         = sces,
                           assay          = assay,
                           slot           = "scale.data",
                           seeds          = 123,
                           ncores         = 4,
                           min.cells      = 3,
                           min.feature    = 0,
                           custom         = T,
                           geneset        = signature_for_deg,
                           msigdb         = T,
                           species        = "Homo sapiens",
                           category       = "H",
                           subcategory    = NULL,
                           geneid         = "symbol",
                           method         = sig_methods,
                           aucell.MaxRank = 2000,
                           ucell.MaxRank  = 2000,
                           kcdf           = 'Gaussian')
      ##########################################

      # 返回一个Seurat对象，富集分数矩阵存放在RNA外的assay中
      Seurat::Assays(sces)
      #> [1] "RNA"       "AUCell"    "UCell"     "singscore" "ssgsea"
      sces@assays$ssgsea
      ##########################################

      # 5.整合差异基因集
      # Wlicox test is perform to all enrichment score matrixes and gene sets
      # with adjusted p value &lt; 0.05 are used to integrated through RRA.
      # Among them, Gene sets with p value &lt; 0.05 are statistically
      # significant and common differential in all gene sets enrichment analysis
      # methods. All results are saved in a list.
      result.dge<-list(NULL)
      # groups<-c("Model1_merge_subcluster", "integrated_snn_res.1")
      groups<-groups[groups%in%colnames(sces@meta.data)]
      names_group<-NULL
      for (j in 1:length(groups)) {

        group<-groups[j]
        #去掉细胞数量很少的cluster==大于50个
        #########################################
        input<-unassign_cell(sce                =  sces,
                             cluster            = group,
                             ignore_cell_prefix = ignore,
                             min_cell_count     = min_cell_count,
                             new_col            = NULL,
                             delete_unassigned  = T,
                             return_meta_data   = T)

        #' 如果目标变量没有亚组直接跳过此次循环====================

        print(table(input[, group]))

        if(length(unique(input[, group]))==1) next

        names_group<-c(names_group, group)
        result.dge[[j]] <- irGSEA.integrate(object   = sces,
                                            group.by = group,
                                            metadata = NULL,
                                            col.name = NULL,
                                            method   = sig_methods)
      }
      names(result.dge)<-names_group
      class(result.dge)

      # 储存数据
      save(sces, result.dge, file = paste0(path_res$abspath, "0-DE-signatues-collection-of-",celltype,".RData"))
      ##########################################

      #'数据可视化==============================

      ##########################################
      # sces<-subset(sces, scpred_seurat_merge!="unassigned")
      # groups<-c("Model1_merge_subcluster", "integrated_snn_res.1")
      # assays<-sig_methods #"ssgsea"



      for (jj in 1:length(sig_methods)) {

        DefaultAssay(sces) <- sig_methods[jj]

        for(ii in 1:length(groups)){
          group<-groups[ii]

          #去掉细胞数量很少的cluster==大于50个
          #########################################
          sces2<-unassign_cell(sce                =  sces,
                               cluster            = group,
                               ignore_cell_prefix = ignore,
                               min_cell_count     = min_cell_count,
                               new_col            = NULL,
                               delete_unassigned  = T,
                               return_meta_data   = FALSE)

          #' 如果目标变量没有亚组直接跳过此次循环====================

          print(table(sces2@meta.data[, group]))
          if(length(unique(sces2@meta.data[, group]))==1) next

          path2<-creat_folder(path_res$folder_name, paste0(ii,"-",group), paste0(jj,"-",sig_methods[jj]))

          # DefaultAssay(sces) <- sig_methods[jj]
          dong_find_markers(sce                      = sces2,
                            assay                    = sig_methods[jj],
                            slot                     = "scale.data",
                            group                    = group,
                            verbose                  = T,
                            feature_type             = "signature",
                            fig.type                 = "pdf",
                            pt.size                  = 0.5,
                            cols                     = cols,
                            palette                  = palette,
                            seed                     = 1234,
                            show_col                 = F,
                            show_genes               = 7,
                            show_genes_pheatmap      = show_feas_pheatmap,
                            hwidth                   = 19,
                            hheight                  = NULL,
                            show_plot                = T,
                            path                     = path2$folder_name,
                            character_limit          = character_limit_heatmap,
                            recluster                = recluster,
                            dims_for_recluster       = 15,
                            resolution_for_recluster = 0.3,
                            assay_for_recluster      =  assay)

        }

      }


    }


  }


}








