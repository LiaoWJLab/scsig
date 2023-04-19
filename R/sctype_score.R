


# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# sctype_score: calculate ScType scores and assign cell types
#
# @params: scRNAseqData - input scRNA-seq matrix (rownames - genes, column names - cells),
# @params: scale - indicates whether the matrix is scaled (TRUE by default)
# @params: gs_sctyper - list of gene sets positively expressed in the cell type
# @params: gs_sctyper2 - list of gene sets that should not be expressed in the cell type (NULL if not applicable)

#' Title
#'
#' @param scRNAseqData
#' @param scaled
#' @param gs_sctyper
#' @param gs_sctyper2
#' @param gene_names_to_uppercase
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sctype_score <- function(scRNAseqData, scaled = !0, gs_sctyper, gs_sctyper2 = NULL, gene_names_to_uppercase = !0, ...){

  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }

  # marker sensitivity
  marker_stat = sort(table(unlist(gs_sctyper)), decreasing = T);
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs_sctyper),1)),
                                  gene_ = names(marker_stat), strings_sctyperAsFactors = !1)

  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }

  # subselect genes only found in data
  names_gs_sctyper_cp = names(gs_sctyper)
  names_gs_sctyper_2_cp = names(gs_sctyper2)
  gs_sctyper = lapply(1:length(gs_sctyper), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs_sctyper[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs_sctyper2 = lapply(1:length(gs_sctyper2), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs_sctyper2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs_sctyper) = names_gs_sctyper_cp
  names(gs_sctyper2) = names_gs_sctyper_2_cp
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs_sctyper)),]

  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
  }

  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs_sctyper),unlist(gs_sctyper2))), ]

  # combine scores
  es = do.call("rbind", lapply(names(gs_sctyper), function(gs_sctypers_){
    sapply(1:ncol(Z), function(j) {
      gs_sctyper_z = Z[gs_sctyper[[gs_sctypers_]], j]; gz_2 = Z[gs_sctyper2[[gs_sctypers_]], j] * -1
      sum_t1 = (sum(gs_sctyper_z) / sqrt(length(gs_sctyper_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) = list(names(gs_sctyper), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

  es.max
}




# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# auto_detect_tissue_type: automatically detect a tissue type of the dataset
#
# @params: path_to_db_file - DB file with cell types
# @params: seuratObject - The Seurat Object from wich to extract the input scRNA-seq matrix (rownames - genes, column names - cells),
# @params: scale - indicates whether the matrix is scaled (TRUE by default)
# @params: assay - e.g. RNA, SCT, integrated

auto_detect_tissue_type <- function(path_to_db_file, seuratObject, scaled, assay = "RNA", ...){

  # get all tissue types in DB
  db_read = openxlsx::read.xlsx(path_to_db_file)
  tissues_ = unique(db_read$tissueType)
  result_ = c()

  for(tissue in tissues_){ print(paste0("Checking...", tissue));

    # prepare gene sets
    gs_sctyper_list = gene_sets_prepare(path_to_db_file, tissue);

    # prepare obj
    if(scaled){
      obj = as.matrix(seuratObject[[assay]]@scale.data)
    } else {
      obj = as.matrix(seuratObject[[assay]]@counts)
    }

    es.max = sctype_score(scRNAseqData = obj, scaled = scaled, gs_sctyper = gs_sctyper_list$gs_sctyper_positive, gs_sctyper2 = gs_sctyper_list$gs_sctyper_negative,
                          marker_sensitivity = gs_sctyper_list$marker_sensitivity, verbose=!0);

    cL_resutls = do.call("rbind", lapply(unique(seuratObject@meta.data$seurat_clusters), function(cl){

      es.max.cl = sort(rowSums(es.max[ ,rownames(seuratObject@meta.data[seuratObject@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
      head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
    }))

    dt_out = cL_resutls %>% group_by(cluster) %>% top_n(n = 1)

    # return mean score for tissue
    result_ = rbind(result_, data.frame(tissue = tissue, score = mean(dt_out$scores)))
  }

  # order by mean score
  result_ = result_[order(-result_$score),]

  # plot
  barplot(height=result_$score, names=result_$tissue, col=rgb(0.8,0.1,0.1,0.6),
          xlab="Tissue", ylab="Summary score",
          main="The higher summary score, the more likely tissue type is")

  result_
}


################################################


# GNU General Public License v3.0 (https://github.com/IanevskiAleksandr/sc-type/blob/master/LICENSE)
# Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#
# Functions on this page:
# gene_sets_prepare: prepare gene sets and calculate marker sensitivity from input Cell Type excel file
#
# @params: path_to_db_file - DB file with cell types
# @cell_type - cell type (e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain)
#

#' gene set prepare
#' Written by Aleksandr Ianevski <aleksandr.ianevski@helsinki.fi>, June 2021
#' Modified by IOBR organization <iobr2019@163.com>, Oct 2022
#'
#' @param gs_sctyper user can provide gene set file manually
#' @param path_to_db_file - DB file with cell types
#' @param cell_type default is null
#' @param tissue_type default is null
#' @param cell_subset default is null
#'
#' @return
#' @export
#'
#' @examples
gene_sets_prepare <- function(gs_sctyper = NULL, path_to_db_file = NULL, cell_type = NULL, tissue_type = NULL, cell_subset = NULL, study = NULL){


  if(!is.null(gs_sctyper)){

    cell_markers<- gs_sctyper

    message(">>>- The input format must be: ")
    exam<- data.frame("tissueType" = "Immune", "cellName" = "T cell", "geneSymbolmore1" = "CD8,CD4", "geneSymbolmore2" = "",
                      "shortName" = "Tcell", "species" = "human/mouse", "studyIndex" = 1, "PublishedYear" = 2022, "tissue" = "stomach")
    exam<-as.tibble(exam)
    print(exam)
    message(">>>- Please double check... ")

  }else{

    if(!is.null(cell_type)){

      message(">>>-- There are 7 options: base, epithelial, myeloid, tcell, bcell, fibroblast, endothelial")
      if(cell_type=="base") sheet =1
      if(cell_type=="epithelial") sheet =2
      if(cell_type=="myeloid") sheet =3
      if(cell_type=="tcell") sheet =4
      if(cell_type=="bcell") sheet =5
      if(cell_type=="fibroblast") sheet =6
      if(cell_type=="endothelial") sheet =7
    }else{
      sheet = 1
      message(">>>-- Default gene sets was selected to estimation: cell_type = 'base'")
    }

    cell_markers = readxl::read_excel(path_to_db_file, sheet = sheet)
  }


  if(!is.null(tissue_type)){

    message(">>>-- Options for tissue_type: ")
    tissue_freq<-as.data.frame(table(cell_markers$tissueType))
    print(tissue_type)


    if(cell_type=="base"){
      cell_markers = cell_markers[cell_markers$tissueType %in%c(tissue_type),]
    }else{
      cell_markers = cell_markers[cell_markers$tissueType %in%c(tissue_type),]
    }


  }


  if(!is.null(cell_subset)){

    message(">>>-- Options for cell_subset: ")
    cell_subset<-as.data.frame(table(cell_markers$subset))
    print(cell_subset)

    cell_markers = cell_markers[cell_markers$subset == cell_subset,]
  }


  if(!is.null(study)){

    message(">>>-- Options for cell_subset: ")
    cell_subset<-as.data.frame(table(cell_markers$studyIndex))
    print(cell_subset)

    cell_markers = cell_markers[cell_markers$studyIndex == study,]
  }

  cell_markers$geneSymbolmore1 = gs_sctyperub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gs_sctyperub(" ","",cell_markers$geneSymbolmore2)

  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){

    markers_all = gs_sctyperub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)

    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })

  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){

    markers_all = gs_sctyperub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)

    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })

  cell_markers$geneSymbolmore1 = gs_sctyperub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gs_sctyperub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gs_sctyperub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gs_sctyperub(" ","",cell_markers$geneSymbolmore2)

  gs_sctyper = lapply(1:nrow(cell_markers), function(j) gs_sctyperub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs_sctyper) = cell_markers$cellName
  gs_sctyper2 = lapply(1:nrow(cell_markers), function(j) gs_sctyperub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs_sctyper2) = cell_markers$cellName

  list(gs_sctyper_positive = gs_sctyper, gs_sctyper_negative = gs_sctyper2)
}




#' Get gene sets from DB, or modified gene sets data
#'
#' @param gs_sctyper user can provide gene set file manually
#' @param path_to_db_file - DB file with cell types
#' @param cell_type default is null
#' @param tissue_type default is null
#' @param cell_subset default is null
#'
#' @return
#' @export
#'
#' @examples
get_gene_sets <- function(gs_sctyper = NULL, path_to_db_file = NULL, cell_type = NULL, tissue_type = NULL, cell_subset = NULL){


  if(!is.null(gs_sctyper)){

    cell_markers<- gs_sctyper

    message(">>>- The input format must be: ")
    exam<- data.frame("tissueType" = "Immune", "cellName" = "T cell", "geneSymbolmore1" = "CD8,CD4", "geneSymbolmore2" = "",
                      "shortName" = "Tcell", "species" = "human/mouse", "studyIndex" = 1, "PublishedYear" = 2022, "tissue" = "stomach")
    exam<-as.tibble(exam)
    print(exam)
    message(">>>- Please double check... ")

  }else{

    if(!is.null(cell_type)){

      message(">>>-- There are 7 options: base, epithelial, myeloid, tcell, bcell, fibroblast, endothelial")
      if(cell_type=="base") sheet =1
      if(cell_type=="epithelial") sheet =2
      if(cell_type=="myeloid") sheet =3
      if(cell_type=="tcell") sheet =4
      if(cell_type=="bcell") sheet =5
      if(cell_type=="fibroblast") sheet =6
      if(cell_type=="endothelial") sheet =7
    }else{
      sheet = 1
      message(">>>-- Default gene sets was selected to estimation: cell_type = 'base'")
    }

    cell_markers = readxl::read_excel(path_to_db_file, sheet = sheet)
  }


  if(!is.null(tissue_type)){

    message(">>>-- Options for tissue_type: ")
    tissue_freq<-as.data.frame(table(cell_markers$tissueType))
    print(tissue_type)

    cell_markers = cell_markers[cell_markers$tissueType == tissue_type,]
  }


  if(!is.null(cell_subset)){

    message(">>>-- Options for cell_subset: ")
    cell_subset<-as.data.frame(table(cell_markers$subset))
    print(cell_subset)

    cell_markers = cell_markers[cell_markers$subset == cell_subset,]
  }


  cell_markers$geneSymbolmore1 = gs_sctyperub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gs_sctyperub(" ","",cell_markers$geneSymbolmore2)

  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){

    markers_all = gs_sctyperub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)

    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })

  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){

    markers_all = gs_sctyperub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)

    if(length(markers_all) > 0){
      markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })

  cell_markers$geneSymbolmore1 = gs_sctyperub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gs_sctyperub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gs_sctyperub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gs_sctyperub(" ","",cell_markers$geneSymbolmore2)

 return(cell_markers)
}

