#' @title Annotate NetBox graph
#'
#' @description
#' This function annotates the graph based on user input. If a table of color codes 
#' for interaction types is provided, then the edges will be colored accordingly by 
#' interaction types. If \code{"directed"} is TRUE, then the edges will be arrows 
#' with the same directionality as the original input network for NetBox. If 
#' \code{"linker"} is TRUE, then linker nodes will be shown as squares while 
#' non-linker nodes stay as circles. 
#'
#' @param netboxResults Output from geneConnector(). 
#'   A list with six lists (i.e. netboxGraph, netboxCommunity, netboxOutput, 
#'                          nodeType, moduleMembership, neighborData)
#'     netboxGraph is an igraph object.
#'     netboxCommunity is an igraph object.
#'     netboxOutput is a data frame.
#'     nodeType is a data frame.
#'     moduleMembership is a data frame.
#'     neighborData is a data frame.
#'
#' @param edgeColors A table (no header) containing hex color codes for interaction types. 
#'                   The first column is interaction type and the second column is hex color code.
#' 
#' @param directed TRUE or FALSE
#'
#' @param linker TRUE or FALSE
#'
#' @return annotated version of netboxGraph
#'
#' @author Guanlan Dong, \email{guanlan_dong@g.harvard.edu}
#'
#' @examples
#' data(netbox2010)
#'
#' sifNetwork<-netbox2010$network
#' graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
#' 
#' geneList<-as.character(netbox2010$geneList)
#' 
#' results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
#'                       pValueAdj='BH',pValueCutoff=0.05,
#'                       communityMethod='lec',keepIsolatedNodes=FALSE)
#' 
#' netboxGraphAnnotated <- annotateGraph(netboxResults = results,
#'                                       edgeColors = "interaction_type.color.txt",
#'                                       directed = TRUE,
#'                                       linker = TRUE)
#' 
#' @concept netboxr
#' @export
#' @import igraph
annotateGraph <- function(netboxResults, edgeColors = NULL, directed = TRUE, linker = TRUE){
  # Extract original edges with interaction types and directions from netboxOutput
  edges <- netboxResults$netboxOutput
  # Reorder columns 
  edges <- edges[, c(1,3,2)]
  # If color codes for interaction types are provided
  if (!is.null(edgeColors)){
    edges$interactionColor <- edgeColors[,2][match(edges$interaction, edgeColors[,1])]
  }
  # Extract nodes from netboxGraph so that they are in the same order to retain correct node colors
  nodes <- data.frame(gene = as_ids(V(netboxResults$netboxGraph)))
  # Add back node type
  nodes$type <- netboxResults$nodeType$type[match(nodes$gene, netboxResults$nodeType$name)]
  # If linker is set to TRUE
  if (linker){
    nodes$shape <- ifelse(nodes$type == "linker", "square", "circle")
  }
  # Re-create the netbox graph
  netboxGraphAnnotated <- graph_from_data_frame(d = edges, vertices = nodes, directed = directed)
  
  return(netboxGraphAnnotated)
}
