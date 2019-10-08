#' @title Simplify sif network into igraph network graph object
#'
#' @description
#' This function removes duplicated edges and loops to create an igraph graph
#' object from tab delimited sif formatted network file.
#'
#' @details
#' For undirected graph, \code{networkSimplify} removes duplicated edges
#' and loops to create an igraph graph object from tab delimited sif
#' formatted network file.
#'
#' For directed graph, \code{networkSimplify} selects the first edge and
#' removes the rest duplicated edges and loops to create an igraph graph
#' object from tab delimited sif
#' formatted network file.
#'
#' @param sifNetwork A file with sif network format (There are three columns
#'        in the file separated by tab, nodeA interactionType nodeB )
#'
#' @param directed Logical, treat network as directed or undirected graph
#'
#' @return a igraph graph object
#'
#' @author Eric Minwei Liu, \email{emliu.research@gmail.com}
#'
#' @examples
#' data(netbox2010)
#'
#' sifNetwork <- netbox2010$network
#' graphReduced <- networkSimplify(sifNetwork, directed = FALSE)
#' @concept netboxr
#' @export
#' @import igraph
networkSimplify <- function(sifNetwork, directed = FALSE) {
  relationsTable <- sifNetwork
  # relationsTable<-read.table('PCWithLoc_PPI_05092014_simple.sif',header=FALSE,sep='\t')
  removeColumn <- c(1, 3)
  relationsTableAttr <- relationsTable[, -removeColumn, drop = FALSE]
  interactions <- data.frame(relationsTable[, 1], relationsTable[, 3], relationsTableAttr, stringsAsFactors = FALSE)

  # load graph as un-directed or directed graph
  graphFull <- graph.data.frame(interactions, directed = directed)
  numOfNodes <- length(V(graphFull))
  numOfEdges <- length(E(graphFull))
  cat(sprintf("Loading network of %s nodes and %s interactions\n", numOfNodes, numOfEdges))
  if (directed == TRUE) {
    directionality <- "directed"
  } else {
    directionality <- "undirected"
  }
  cat(sprintf("Treated as %s network \n", directionality))

  # print(graphFull, e=TRUE, v=TRUE) The opposite operation get.data.frame(graphFull, what='vertices')
  # get.data.frame(graphFull, what='edges')

  # remove multiple interactions among the same pair of nodes and remove interaction loops
  graphReduced <- simplify(graphFull,
    remove.multiple = TRUE, remove.loops = TRUE,
    edge.attr.comb = "first"
  )
  numOfNodes <- length(V(graphReduced))
  numOfEdges <- length(E(graphReduced))
  cat(sprintf("Removing multiple interactions and loops\n"))
  cat(sprintf("Returning network of %s nodes and %s interactions\n", numOfNodes, numOfEdges))

  return(graphReduced)
}
