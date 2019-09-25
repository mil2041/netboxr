#' network coming with PLOS ONE paper.
#'
#' Loading PC2V7_03052015_woHPRD containing 19147 nodes and 2267149 interactions.
#' Treated as undirected network. 
#' Removing multiple interactions and loops.
#' Returning igraph network of 19147 nodes and 2196039 interactions.
#'
#' @format A data frame with 19147 nodes and 2196039 interactions:
#' \describe{
#'   \item{name}{vertex gene name}
#'   \item{edges}{interaction types}
#'   ...
#' }
#' @source \url{http://www.pathwaycommons.org/pc2/downloads}
#' 
#' @return a data.frame
#' 
#' @concept netboxr
"netbox2010"

