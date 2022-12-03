



##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  pkgVersion <- packageDescription(pkgname, fields="Version")
  msg <- paste0("==========================================================================\n",
                " ", pkgname, " v", pkgVersion, "  ",

                "  For help: https://github.com/DongqiangZeng0808/", pkgname, "\n\n")

  citation <-paste0(" If you use ", pkgname, " in published research, please cite:\n",
                    " Wellcome to my Single cell data processing.\n",

                    " Frontiers in Immunology, 2050, 7(5), 737-750", "\n",
                    " DOI: 10.1158/2326-6066.nejm-18-2050 ","\n" ,
                    " PMID: 3838438","\n",
                    "===========================================================================")

  packageStartupMessage(paste0(msg, citation))
}




# .onLoad <- function(libname, pkgname) {
#   op <- options()
#   op.devtools <- list(
#     devtools.path = "~/R-dev",
#     devtools.install.args = "",
#     devtools.name = "DongqiagnZeng0808",
#     devtools.desc.author = '"Dongqiang Zeng <dognqiangzeng0808@gmail.com> [aut, cre]"',
#     devtools.desc.license = "GPL-3",
#     devtools.desc.suggests = NULL,
#     devtools.desc = list()
#   )
#   toset <- !(names(op.devtools) %in% names(op))
#   if(any(toset)) options(op.devtools[toset])
#
#   requireNamespace(c('GSVA','limma','DESeq2','tidyverse','MASS',"ggplot2"))
#
#   invisible()
# }
