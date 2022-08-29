#' @title Annotate NetBox graph
#'
#' @description
#' This function annotates the graph based on user input. 
#' 
#' @details 
#' If a table of color codes for interaction types is provided, then the edges 
#' will be colored accordingly by interaction types. 
#' If \code{directed} is TRUE, then the edges will be arrows 
#' with the same directionality as the original input network for NetBox. If 
#' \code{linker} is TRUE, then linker nodes will be shown as squares while 
#' non-linker nodes stay as circles. 
#'
#' @param netboxResults Output from \code{geneConnector} function. 
#' a list of returned netboxr results 
#' * netboxGraph: igraph object of NetBox algorithm identified network nodes 
#' and connections
#' * netboxCommunity: igraph object of network community assignment
#' * netboxOutput: data frame of NetBox algorithm identified network nodes 
#' and connections
#' * nodeType: data frame of node types ("candidate" or "linker") 
#' in the NetBox algorithm indentified network.
#' * moduleMembership: data frame of module (community) membership.
#' * neighborData: data frame of information of nodes directly connected to 
#' candidate gene nodes.
#' @md
#'
#' @param edgeColors table containing hex color codes for interaction types. 
#'                   The first column is interaction type 
#'                   and the second column is hex color code.
#' 
#' @param directed boolean value indicating whether the NetBox algorithm 
#' identified network is directed or undirected (default = FALSE)
#'
#' @param linker boolean value indicating whether "linker" nodes exist 
#' in the NetBox algorithm identified network or not (default = TRUE)
#'
#' @return annotated version of netboxGraph
#'
#' @author Guanlan Dong, \email{guanlan_dong@g.harvard.edu}
#'
#' @examples
#' data(netbox2010)
#' interaction_type_color <- read.csv(system.file("interaction_type.color.txt", package = "netboxr"),
#'                                    header=TRUE, sep="\t", stringsAsFactors=FALSE)
#' 
#' sifNetwork<-netbox2010$network
#' graphReduced <- networkSimplify(sifNetwork,directed = FALSE)
#' 
#' geneList <- netbox2010$geneList[1:100]
#' 
#' results <- geneConnector(geneList = geneList, networkGraph = graphReduced, 
#'                          directed = FALSE, pValueAdj = "BH", pValueCutoff = 2e-5, 
#'                          communityMethod = "lec", keepIsolatedNodes = FALSE)
#'
#' netboxGraphAnnotated <- annotateGraph(netboxResults = results,
#'                                       edgeColors = interaction_type_color,
#'                                       directed = TRUE,
#'                                       linker = TRUE)
#' 
#' # As an example, plot both the original and the annotated graphs
#' ll <- layout_with_fr(results$netboxGraph) # Save the layout for easier comparison
#' # Plot original graph
#' \dontrun{
#' pdf("originalGraph.pdf", width = 50, height = 50)
#' plot(results$netboxCommunity, results$netboxGraph, layout = ll,
#'      vertex.size=3)
#' dev.off()
#' # Plot annotated graph
#' pdf("annotatedGraph.pdf", width = 50, height = 50)
#' plot(results$netboxCommunity, netboxGraphAnnotated, layout = ll,
#'      vertex.size = 3,
#'      vertex.shape = V(netboxGraphAnnotated)$shape,
#'      edge.color = E(netboxGraphAnnotated)$interactionColor,
#'      edge.width = 3)
#' # Add legend
#' ind <- which(interaction_type_color$INTERACTION_TYPE %in% E(netboxGraphAnnotated)$interaction)
#' legend_interaction_type <- interaction_type_color$INTERACTION_TYPE[ind]
#' legend_interaction_type_color <- interaction_type_color$COLOR[ind]
#' legend(x=-1.1, y=1.1, 
#'        legend=c("Candidate", "Linker"),
#'        pch=c(19, 15), # solid circle, filled square
#'        pt.cex = 8,
#'        bty="n",
#'        title="Node Types",
#'        cex=4, ncol=1)
#' legend(x=-1.15, y=0.95, 
#'        legend=legend_interaction_type,
#'        col = legend_interaction_type_color,
#'        lty = 1, lwd = 10,
#'        bty="n",
#'        title="Interaction Types (Edges)",
#'        cex=4, ncol=1)
#' dev.off()
#' }
#' 
#' @concept netboxr
#' @export
#' @import igraph
annotateGraph <- function(netboxResults, edgeColors = NULL, directed = FALSE, linker = TRUE){
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
