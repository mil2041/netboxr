#' network coming with Cerami et al. PLoS One 2010 paper.
#'
#' Loading netbox2010 containing 9264 nodes and 68111 interactions.
#' Treated as undirected network.
#' After removing multiple interactions and loops.
#' Returning igraph network of 9264 nodes and 68111 interactions.
#'
#' @format A data frame with 9264 nodes and 68111 interactions:
#' \describe{
#'   \item{name}{vertex gene name}
#'   \item{edges}{interaction types}
#'   ...
#' }
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/20169195}
#'
#' @return a data.frame
#'
#' @concept netboxr
"netbox2010"

#' Pathway Commons V8 Reactome
#' 
#' Contains an example gene list and Pathway Commons V8 Reactome dataset for annotateGraph().
#'
#' @format A list of 354 genes and a data frame of 246590 interactions 
#' \describe{
#'   \item{geneList}{an example list of genes}
#'   \item{network}{Pathway Commons V8 Reactome}
#'   ...
#' }
#' 
#' @source \url{https://www.pathwaycommons.org}
#'
#' @return A list of two elements.
#'
#' @concept netboxr
"pathway_commons_v8_reactome"