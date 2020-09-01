#' @title Generate sub-network mapping from a list of candidate genes
#' 
#' @description 
#' This function generates sub-network mapping from a list of candidate genes
#' 
#' @details 
#' P-value correction methods include the Bonferroni correction 
#' ("bonferroni") or Benjamini & Hochberg ("BH"). Community detection methods 
#' include using edge betweeness score ("ebc") or using leading eigenvector 
#' method ("lec) 
#'
#' @param geneList character vector containing a list of candidate genes
#' @param networkGraph igraph network graph object. This igraph object contains
#' curated network information 
#' @param directed boolean value indicating whether the input network is
#' directed or undirected (default = FALSE)
#' @param pValueAdj string for p-value correction method c("BH", "Bonferroni")
#' as described in the details section (default = "BH")
#' @param pValueCutoff numeric value of p-value cutoff for linker nodes 
#' (default = 0.05)
#' @param communityMethod string for community detection method c("ebc","lec")
#' as described in the details section (default = "ebc")
#' @param keepIsolatedNodes A boolean value indicating whether to keep isolated
#' nodes in the netboxr result (default = FALSE)
#'
#' @return a list of returned netboxr results 
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
#' @author Eric Minwei Liu, \email{emliu.research@gmail.com}
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
#'
#' names(results)
#' 
#' plot(results$netboxGraph, layout = layout_with_fr)
#' 
#'
#' write.table(results$netboxOutput,
#'   file = "network.sif", sep = "	",
#'   quote = FALSE, col.names = FALSE, row.names = FALSE
#' )
#'
#' write.table(results$neighborData,
#'   file = "neighborList.txt", sep = "	",
#'   quote = FALSE, col.names = TRUE, row.names = FALSE
#' )
#'
#' write.table(results$moduleMembership,
#'   file = "memb.ebc.txt", sep = "	",
#'   quote = FALSE, col.names = FALSE, row.names = FALSE
#' )
#' #
#' write.table(results$nodeType,
#'   file = "nodeType.txt", sep = "	", quote = FALSE,
#'   col.names = FALSE, row.names = FALSE
#' )
#' #
#' 
#' @concept netboxr
#' 
#' @export
#' @import igraph
#' @import jsonlite
#' @import parallel
#' @import data.table
#' @import gplots
#' @import plyr
#' @importFrom stats p.adjust
#' @importFrom stats phyper
#' @importFrom clusterProfiler enrichGO
#' @importFrom clusterProfiler bitr
#' @importFrom DT datatable
geneConnector <- function(geneList, networkGraph, directed = FALSE,
                          pValueAdj = "BH", pValueCutoff = 0.05,
                          communityMethod = "ebc", keepIsolatedNodes = FALSE) {
  
  
  pValueAdj<-match.arg(pValueAdj)
  
  graphReduced <- networkGraph

  # Get the genes from the input that overlap with the vertexs in the network
  geneOverlap <- NULL
  geneOverlap$name <- intersect(geneList, V(graphReduced)$name)
  geneOverlap$idx <- match(geneOverlap$name, V(graphReduced)$name)
  geneOverlap$type <- rep("candidate", length(geneOverlap$name))
  message(sprintf(
    "%s / %s candidate nodes match the name in the network of %s 
                nodes \n",
    length(geneOverlap$name), length(geneList), length(V(graphReduced)$name)
  ))

  neighborList <- NULL
  neighborList$numOfgraphReducedGene <- length(V(graphReduced))
  neighborList$numOfgeneOverlap <- length(geneOverlap$idx)


  neighborTotal <- NULL
  neighborTemp <- NULL
  numOfgeneOverlap<-length(geneOverlap$idx)
  
  for (k in seq_len(numOfgeneOverlap)) {
    neighborTemp <- neighbors(graphReduced, v = geneOverlap$idx[k], mode = "all")
    neighborTotal <- union(neighborTotal, neighborTemp)
  }
  neighborList$reducedGraphId <- setdiff(neighborTotal, geneOverlap$idx)


  neighborList$name <- NULL
  neighborList$localDegree <- NULL
  neighborList$pValueRaw <- NULL
  neighborList$oddsRatio <- NULL

  numOfreducedGraphID<-length(neighborList$reducedGraphId)
  
  for (i in seq_len(numOfreducedGraphID)) {
    neighborList$name[i] <- get.vertex.attribute(graphReduced,
      name = "name",
      index = neighborList$reducedGraphId[i]
    )
  }

  neighborList$globalDegree <- degree(graphReduced, v = neighborList$reducedGraphId, mode = "all")

  for (j in seq_len(numOfreducedGraphID)) {
    globalNeighborId <- neighbors(graphReduced, v = neighborList$reducedGraphId[j], mode = 1)

    localNeighborId <- intersect(geneOverlap$id, globalNeighborId)

    neighborList$localDegree[j] <- length(localNeighborId)

    # a: number of nodes A node links to the gene list in the graph
    # b: number of nodes A node links in the graph globally
    # c: total number of nodes in graph minus total number of nodes that A node links globally
    # d: total number of nodes in the gene list in the graph

    a <- neighborList$localDegree[j]
    b <- neighborList$globalDegree[j]
    c <- (neighborList$numOfgraphReducedGene - neighborList$globalDegree[j])
    d <- neighborList$numOfgeneOverlap

    # hypergeometric distribution for p-value calculation neighborList$pValueRaw[j]<-(1-phyper((a-1),b,c,d))

    # http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
    # https://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html

    neighborList$pValueRaw[j] <- phyper((a - 1), b, c, d, lower.tail = FALSE)

    # https://en.wikipedia.org/wiki/Odds_ratio
    # https://www.researchgate.net/post/How_to_calculate_OR_odd_ratio_if_one_of_groups_is_0_in_a_case-control_study
    # Add 0.5 to each cell to avoid zero value cell problem
    correctionTerm <- 0.5

    n11 <- neighborList$localDegree[j] + correctionTerm
    n10 <- neighborList$globalDegree[j] - neighborList$localDegree[j] + correctionTerm
    n01 <- neighborList$numOfgeneOverlap - neighborList$localDegree[j] + correctionTerm
    n00 <- neighborList$numOfgraphReducedGene - neighborList$globalDegree[j] - n01 + correctionTerm

    neighborList$oddsRatio[j] <- (log(n11) + log(n00) - log(n10) - log(n01))
  }

  neighborListFrame <- data.frame(neighborList$reducedGraphId, neighborList$name, neighborList$localDegree,
    neighborList$globalDegree, neighborList$pValueRaw, neighborList$oddsRatio,
    stringsAsFactors = FALSE
  )
  colnames(neighborListFrame) <- c("idx", "name", "localDegree", "globalDegree", "pValueRaw", "oddsRatio")
  neighborListFrame <- neighborListFrame[order(neighborListFrame$pValueRaw), ]

  localDegreeCutoff <- 2
  neighborListFrame <- neighborListFrame[neighborListFrame$localDegree >= localDegreeCutoff, ]

  message(sprintf("Only test neighbor nodes with local degree equals or exceeds %s\n", localDegreeCutoff))
  message(sprintf("Multiple hypothesis corrections for %s neighbor nodes in the network\n", nrow(neighborListFrame)))

  neighborListFrame$pValueFDR <- p.adjust(neighborListFrame$pValueRaw, method = pValueAdj)
  
  linkerListFrame <- neighborListFrame[neighborListFrame$pValueFDR < pValueCutoff, ]
 
  message(sprintf(
    "For p-value %s cut-off, %s nodes were included as linker nodes\n",
    pValueCutoff, dim(linkerListFrame)[1]
  ))
  linkerListFrame$type <- rep("linker", dim(linkerListFrame)[1])

  # include linker nodes into the list of initial candidate nodes
  selectedGene <- NULL
  selectedGene$idx <- union(geneOverlap$idx, linkerListFrame$idx)
  selectedGene$type <- c(geneOverlap$type, linkerListFrame$type)
  selectedGene$name <- c(geneOverlap$name, as.character(linkerListFrame$name))
  selectedNodeType <- data.frame(selectedGene$name, selectedGene$type, stringsAsFactors = FALSE)
  colnames(selectedNodeType) <- c("name", "type")



  message(sprintf(
    "Connecting %s candidate nodes and %s linker nodes\n",
    length(geneOverlap$name), dim(linkerListFrame)[1]
  ))
  # sub-network of linker nodes and candidate nodes
  graphTemp <- induced.subgraph(graphReduced, selectedGene$idx)

  # add vertex attributes
  graphTemp <- set_vertex_attr(graphTemp, name = "nodeType", index = V(graphTemp), value = "candidate")
  graphTemp <- set_vertex_attr(graphTemp,
    name = "nodeType", index = as.character(linkerListFrame$name),
    value = "linker"
  )


  # remove isolated vertex
  isolatedNodes <- names(which(degree(graphTemp) < 1))

  if (keepIsolatedNodes) {
    graphOutput <- graphTemp
    message(sprintf("Keep %s isolated candidate nodes from the input\n", length(isolatedNodes)))
  } else {
    graphOutput <- delete.vertices(graphTemp, which(degree(graphTemp) < 1))
    message(sprintf("Remove %s isolated candidate nodes from the input\n", length(isolatedNodes)))
  }


  numOfNodes <- length(V(graphOutput))
  numOfEdges <- length(E(graphOutput))

  message(sprintf("Final network contains %s nodes and %s interactions\n", numOfNodes, numOfEdges))

  # Assign community membership for the sub-network

  if (communityMethod == "ebc") {
    message(sprintf("Detecting modules using \"edge betweeness\" method\n"))
    community <- edge.betweenness.community(graphOutput)
    moduleMembership <- membership(community)
  }

  if (communityMethod == "lec") {
    message(sprintf("Detecting modules using \"leading eigenvector\" method\n"))
    community <- leading.eigenvector.community(graphOutput, options = list(maxiter = 1e+06))
    moduleMembership <- membership(community)
  }

  moduleMembership <- moduleMembership[names(moduleMembership)]
  moduleMembershipFrame <- data.frame(names(moduleMembership), moduleMembership, stringsAsFactors = FALSE)
  colnames(moduleMembershipFrame) <- c("geneSymbol", "membership")
  moduleMembershipFrame <- moduleMembershipFrame[order(moduleMembershipFrame$membership), ]


  netbox <- NULL
  
  netbox$edgelist <- get.edgelist(graphOutput)

  edgeLabelList <- get.edge.attribute(graphOutput)$V2

  # Added for Pathway Commons datasets
  if (is.null(edgeLabelList)) {
    edgeLabelList <- get.edge.attribute(graphOutput)$INTERACTION_TYPE
  }

  mergedEdgeLabels <- lapply(edgeLabelList, function(x) {
    content <- unlist(unique(strsplit(x, ";")))
    contentMerged <- paste(content, collapse = ";")
    numOfcontent <- length(content)
    data.frame(numOfcontent, contentMerged, stringsAsFactors = FALSE)
  })

  mergedEdgeLabels <- do.call(rbind, mergedEdgeLabels)

  netboxOutput <- data.frame(netbox$edgelist[, 1], mergedEdgeLabels$contentMerged, netbox$edgelist[, 2],
    stringsAsFactors = FALSE
  )
  colnames(netboxOutput) <- c("geneA", "interaction", "geneB")

  if (keepIsolatedNodes) {
    netboxOutputExtra <- data.frame(isolatedNodes, rep("INTERACT", length(isolatedNodes)), isolatedNodes,
      stringsAsFactors = FALSE
    )
    colnames(netboxOutputExtra) <- c("geneA", "interaction", "geneB")
    netboxOutput <- rbind(netboxOutput, netboxOutputExtra)
  }

  result <- list(
    netboxGraph = graphOutput, 
    netboxCommunity = community,
    netboxOutput = netboxOutput, 
    nodeType = selectedNodeType,
    moduleMembership = moduleMembershipFrame, 
    neighborData = neighborListFrame
  )

  return(result)
}
