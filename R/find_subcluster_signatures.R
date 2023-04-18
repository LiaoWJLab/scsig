




<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
#' Gene-set expression calculation and visualization for cell subtypes
#'
#' For the cell type of interest, firstly, `unassign_cell()` is used to remove the cell subtypes with few cells.
#' Secondly, use `irGSEA.score()` to calculate gene-set enrichment scores score by optional methods
#' and `irGSEA.integrate` to integrate differential gene set calculated by all enrichment score matrixes.
#' See more details in [irGSEA.score](http://127.0.0.1:56860/help/library/scsig/help/irGSEA.score) and [irGSEA.integrate](http://127.0.0.1:56860/help/library/scsig/help/irGSEA.integrate.).
#' Finally, data visualization is achieved by `dong_find_markers()`.
#'
#' @param sce Seurat object
#' @param assay Assay to pull from, e.g. RNA, SCT, integrated
#' @param col_celltype Name of column where cell types in the metadata
#' @param col_sub_celltype  Name of column where cell subtypes in the metadata
#' @param groups Name of one or more metadata columns to group cells by
#' @param signature_for_deg Gene-set list
#' @param min_cell_count Minimal cell counts in each cluster, default is 50
#' @param no_limitation_celltypes A character vector corresponding to cell types with no cell count limitation
#' @param sig_methods Method used to estimate Gene-set expression score, default is `PCAscore`,`ssgsea`, `AUCell`
#' @param path Path of the output saving directory
#' @param sig_file_name Name of the storage file
#' @param index User can choose specific cell types to save computing source, default is NULL
#' @param cols Vector of colors, users can define the cols manually. This may also be a single character, such as normal and random,  to a palette as specified by `palettes(){IOBR}`
#'             See [palettes](http://127.0.0.1:60491/help/library/IOBR/html/palettes.html) for details
#' @param show_feas_pheatmap Number of genes displayed in the heatmap. Default is 8
#' @param character_limit_heatmap Limit length of character. Default is 50
#' @param remove_other_celltypes Default is TURE
#' @param find_cluster_sig If TRUE, the analysis process will be performed in all clusters. Default is FALSE
#' @param groups_for_cluster Name of one or more metadata columns to group cells by. The parameter works when find_cluster_sig = TRUE
#' @param signature_character_limit Limit length of gene-set name, default is 60
#' @param seed Seed of the random number generator, default is 123. The parameter works when cols ="random"
#' @param palette Numeric value corresponding with color palette. Default is 1, other options: 2, 3, 4
#' @param show_col Whether to show color palettes
#'
#' @return Seurat Object
#' @export
#' @author Dongqiang Zeng
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
  # #select groups
  # groups<-c("Model1_merge_subcluster", "SCT_snn_res.1")
  # signature_for_deg<-signature_collection
  # sig_file_name<- "mysig"
  # sig_methods<-c("PCAscore","ssgsea", "AUCell")

  includ_sig<-nchar(names(signature_for_deg)) <= signature_character_limit
  signature_for_deg<- signature_for_deg[includ_sig]

  #creat path to save results
  file_name<-creat_folder(paste0(path, sig_file_name))


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

  #identify differential signatures of different cell types#################################################################################


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

    # return a Seurat object including score matrix.
    Seurat::Assays(sce_a)
    #> [1] "RNA"       "AUCell"    "UCell"     "singscore" "ssgsea"
    # sce_a@assays$ssgsea
    ##########################################

    # integrat DEG
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
      #Remove the cluster with less than 50 cells
      #########################################
      input<-unassign_cell(sce                = sce_a,
                           cluster            = group,
                           ignore_cell_prefix = NULL,
                           min_cell_count     = min_cell_count,
                           new_col            = NULL,
                           delete_unassigned  = T,
                           return_meta_data   = T)

      #' If the target variable has no subgroup, skip the loop directly====================

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

    # save data
    save(sce_a, result.dge, file = paste0(path_res$abspath, "0-DE-signatues-of-","all-clusters",".RData"))
    ##########################################


    for (jj in 1:length(sig_methods)) {

      DefaultAssay(sce_a) <- sig_methods[jj]

      for(ii in 1:length(groups_for_cluster)){
        group<-groups_for_cluster[ii]

        #Remove the clusters with less than 50 cells
        #########################################
        sces2<-unassign_cell(sce                = sce_a,
                             cluster            = group,
                             ignore_cell_prefix = NULL,
                             min_cell_count     = min_cell_count,
                             new_col            = NULL,
                             delete_unassigned  = T,
                             return_meta_data   = FALSE)

        #' If the target variable has no subgroup, skip the loop directly====================

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
        dong_find_markers(sce        = sces2,
                          assay      = sig_methods[jj],
                          slot       = "scale.data",

                          group      = group,
                          verbose    = T,
                          feature_type = "signature",
                          fig.type   = "pdf",
                          pt.size    = 0.5,
                          cols       = cols,
                          seed       = seed,
                          palette    = palette,
                          show_col   = show_col,
                          show_genes = 7,
                          show_genes_pheatmap = show_feas_pheatmap,
                          hwidth     = 19,
                          hheight    = NULL,
                          show_plot  = T,
                          path       = path2$folder_name,
                          character_limit = character_limit_heatmap,
                          recluster = recluster,
                          dims_for_recluster = 8,
                          resolution_for_recluster = 0.3,
                          assay_for_recluster =  assay)

      }

    }

  }


  ###############################++++find_subcluster_signatures####################
  #identify differential signatures of different cell subtypes
  if(!is.null(col_sub_celltype)){

    for(i in index){

      celltype<-celltypes[i]
      message(paste0(">>>>----Processing celltypes ", celltype))
      ############################################


      if(!is.null(assay)) DefaultAssay(sce)<-assay

      Idents(sce)<- col_celltype

      sces<-subset(sce, idents = celltype)


      if(remove_other_celltypes){
        sces@meta.data[, col_sub_celltype]<-ifelse(grepl(sces@meta.data[, col_sub_celltype],pattern = celltype), sces@meta.data[, col_sub_celltype], "unassigned" )
      }

      # sces<-subset(sce, Model1_merge_no_rejection == celltype)
      ###########################################

      #Remove the clusters with less than 50 cells
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

      #'  If the target variable has no subgroup, skip the loop directly====================
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
      # calculate enrichment score
      # When ncore setting is more than 1, the following error occurs: Error (Valid 'mctype': 'snow' or 'doMC'),
      # You should check your AUCell version to ensure that the version is greater than or equal to 1.14.
      # If you are lazy, you can set ncore to 1 directly, but the running speed will be a little slower.

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

      # Returns a Seurat object, and the enrichment fraction matrix is stored in the assay outside RNA
      Seurat::Assays(sces)
      #> [1] "RNA"       "AUCell"    "UCell"     "singscore" "ssgsea"
      sces@assays$ssgsea
      ##########################################

      # integrat DEG
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
        #Remove the clusters with less than 50 cells
        #########################################
        input<-unassign_cell(sce                =  sces,
                             cluster            = group,
                             ignore_cell_prefix = ignore,
                             min_cell_count     = min_cell_count,
                             new_col            = NULL,
                             delete_unassigned  = T,
                             return_meta_data   = T)

        #'  If the target variable has no subgroup, skip the loop directly====================

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

      # save data
      save(sces, result.dge, file = paste0(path_res$abspath, "0-DE-signatues-collection-of-",celltype,".RData"))
      ##########################################

      #'data visualization==============================

      ##########################################
      # sces<-subset(sces, scpred_seurat_merge!="unassigned")
      # groups<-c("Model1_merge_subcluster", "integrated_snn_res.1")
      # assays<-sig_methods #"ssgsea"



      for (jj in 1:length(sig_methods)) {

        DefaultAssay(sces) <- sig_methods[jj]

        for(ii in 1:length(groups)){
          group<-groups[ii]

          #Remove the clusters with less than 50 cells
          #########################################
          sces2<-unassign_cell(sce                =  sces,
                               cluster            = group,
                               ignore_cell_prefix = ignore,
                               min_cell_count     = min_cell_count,
                               new_col            = NULL,
                               delete_unassigned  = T,
                               return_meta_data   = FALSE)

          #' If the target variable has no subgroup, skip the loop directly====================

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








