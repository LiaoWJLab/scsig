




#' Get gene sets from MSigDB
#'
#' Retrieve a list or data frame of gene sets and their member genes.
#' Users can not only download gene sets from MSigSB online directly,
#' but also get them from local datasets by setting the parameter `msigdb`.
#'
#' @param org Species name of organisms, choose from 'hsa' and 'mus'
#' @param category MSigDB collection abbreviation, such as H or C1
#' @param subcategory MSigDB sub-collection abbreviation, such as CGP or BP
#' @param msigdb Local data frame of gent sets. If NULL, it will download from MSigSB automatically
#' @param format Format of the output. Default is list, otherwise data frame
#'
#' @return A list of gene sets or a data frame of gene sets
#' @export
#'
#' @examples
#' ## get human Hallmark gene sets by default
#' gs<- get_msigdb_gs()
#' ## NOTRUN
#' ## get human Hallmark gene sets from local datasets
#' # load local datasets
#' load("msigdb.RData")
#' ##get human C5 (ontology gene sets) GO:BP (GO Biological Process ontology) gene sets
#' term_genes <- get_msigdb_gs(msigdb = msigdb, org = "hsa", category = "C5", subcategory = "GO:BP", format = "list")
get_msigdb_gs<-function(msigdb = NULL, org = "hsa", category = NULL, subcategory = NULL, format = "list"){
  if(org=="hsa"){
    species<- "Homo sapiens"
  }else if(org=="mus"){
    species<- "Mus musculus"
  }
  ################################################

  if(is.null(category)){
    message(">>>---Category is NULL, default is Hallmark gene sets...")
    category = "H"
  }
  ##################################################
  message(">>>---Categories that can be choosed... ")

  if(is.null(msigdb)){
    m_df = msigdbr::msigdbr(species = species)
    # table(m_df$gs_cat, m_df$gs_subcat)
    a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
    print(as.data.frame(a))

    ##################################################
    term2genes <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)

    term2genes <- term2genes %>% dplyr::select(gs_name, entrez_gene, gene_symbol) %>% as.data.frame()
    ##############################################
  }else{
    m_df = msigdb
    # table(m_df$gs_cat, m_df$gs_subcat)
    a <- m_df %>% dplyr::distinct(gs_cat, gs_subcat) %>% dplyr::arrange(gs_cat, gs_subcat)
    print(as.data.frame(a))

    # category<- sym(!category)
    # subcategory<- sym(!subcategory)
    ############################################
    if(is.null(subcategory)){
      term2genes <- msigdb %>%
        filter(gs_cat == rlang::sym(category)) %>%
        # filter(gs_subcat == rlang::sym(subcategory)) %>%
        dplyr::select(gs_name, entrez_gene, gene_symbol) %>%
        as.data.frame()
      ############################################
    }else{
      term2genes <- msigdb %>%
        filter(gs_cat == rlang::sym(category)) %>%
        filter(gs_subcat == rlang::sym(subcategory)) %>%
        dplyr::select(gs_name, entrez_gene, gene_symbol) %>%
        as.data.frame()
      ############################################
    }

  }

  if(format=="list"){
    term2genes <- select(term2genes, gs_name, gene_symbol) %>%
      as.data.frame %>%
      split(., .$gs_name) %>%
      lapply(., function(x)(x$gene_symbol))
  }

  return(term2genes)
}

####################################################
# gs<- get_msigdb_gs(msigdb = msigdb)
# gs
# ##################################################
# help("calculate_sig_score")
# score<-IOBR::calculate_sig_score(eset = eset_stad, signature = gs, method = "ssgsea")
# score
####################################################

