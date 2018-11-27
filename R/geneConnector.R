#' Generate sub-network mapping from a candidate gene list
#' 
#' @param geneList A vector containing candidate gene list 
#' @param networkGraph An igraph graph object
#' @param directed TRUE of FALSE
#' @param pValueAdj A string for p-value correction method c("BH, "Bonferroni")
#' @param pValueCutoff A number for p-value cutoff for linker nodes
#' @param communityMethod A string for community detection method c("ebc","lec")
#' @param keepIsolatedNodes logic value
#' 
#' @return a list with four lists (i.e. netboxOutput, nodeType, moduleMembership, neighborData)
#'   netboxGraph is igraph object.
#'   netboxCommunity is igraph object.
#'   netboxOutput is a data frame.
#'   nodeType is a data frame.
#'   moduleMembership is a data frame.
#'   neighborData is a data frame.
#'   
#' @author Eric Minwei Liu, \email{emliu.research@gmail.com}
#'  
#' @examples
#' data(netbox2010)
#'
#' geneList<-netbox2010$geneList 
#' sifNetwork<-netbox2010$network
#' graphReduced<-networkSimplify(sifNetwork,directed = FALSE)      
#' result<-geneConnector(geneList=geneList,networkGraph=graphReduced,directed=FALSE,pValueAdj="BH",
#'            pValueCutoff=0.05,communityMethod="lec",keepIsolatedNodes=FALSE)
#'
#' names(result)
#'
#' \dontrun{ plot(result$netboxGraph, layout=layout_with_fr) }
#' 
#' write.table(result$netboxOutput,file="network.sif",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
#' write.table(result$neighborData,file="neighborList.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#' write.table(result$moduleMembership,file="memb.ebc.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
#' write.table(result$nodeType,file="nodeType.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
#' 
#' @concept netboxr
#' @export
#' @import igraph
geneConnector<-function(geneList,networkGraph,directed=FALSE,pValueAdj="BH",pValueCutoff=0.05,communityMethod="lec",keepIsolatedNodes=FALSE){
  
  graphReduced<-networkGraph
  
  # Get the genes from the input that overlap with the vertexs in the network  
  geneOverlap<-{}  
  geneOverlap$name<-intersect(geneList,V(graphReduced)$name)
  geneOverlap$idx<-match(geneOverlap$name,V(graphReduced)$name)
  geneOverlap$type<-rep("candidate",length(geneOverlap$name))
  cat(sprintf("%s / %s candidate nodes match the name in the network of %s nodes \n",length(geneOverlap$name),length(geneList),length(V(graphReduced)$name)))
  
  neighborList<-{}
  neighborList$numOfgraphReducedGene<-length(V(graphReduced))
  neighborList$numOfgeneOverlap<-length(geneOverlap$idx)  
  
  #neighborList$reducedGraphId<-neighbors(graphReduced,v=geneOverlap$idx,mode=1)
  
  neighborTotal<-{}
  neighborTemp<-{}
  for (k in 1:length(geneOverlap$idx))
  {
    neighborTemp<-neighbors(graphReduced,v=geneOverlap$idx[k],mode="all")
    neighborTotal<-union(neighborTotal,neighborTemp)
  }
  neighborList$reducedGraphId<-setdiff(neighborTotal, geneOverlap$idx)
  
  
  neighborList$name<-{}
  neighborList$localDegree<-{}
  neighborList$pValueRaw<-{}
  neighborList$oddsRatio<-{}
  
  for (i in 1:length(neighborList$reducedGraphId))
  {
    neighborList$name[i]<-get.vertex.attribute(graphReduced,name="name",index=neighborList$reducedGraphId[i])  
  }  
  
  neighborList$globalDegree<-degree(graphReduced,v=neighborList$reducedGraphId,mode="all")
  
  for (j in 1:length(neighborList$reducedGraphId))
  {  
    
    globalNeighborId<-neighbors(graphReduced,v=neighborList$reducedGraphId[j],mode=1)
    localNeighborId<-intersect(geneOverlap$id,globalNeighborId)
    neighborList$localDegree[j]<-length(localNeighborId)
    
    # a: number of nodes A node links to the gene list in the graph
    # b: number of nodes A node links in the graph globally
    # c: total number of nodes in graph minus toal number of nodes that A node links globally
    # d: total number of nodes in the gene list in the graph
    
    a<-neighborList$localDegree[j]
    b<-neighborList$globalDegree[j]
    c<-(neighborList$numOfgraphReducedGene - neighborList$globalDegree[j])  
    d<-neighborList$numOfgeneOverlap 
    
    # hypergeometric distribution for p-value calculation
    #neighborList$pValueRaw[j]<-(1-phyper((a-1),b,c,d))
    
    # http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
    # https://mengnote.blogspot.com/2012/12/calculate-correct-hypergeometric-p.html
    neighborList$pValueRaw[j]<-phyper((a-1),b,c,d, lower.tail = FALSE)
    
    # https://en.wikipedia.org/wiki/Odds_ratio
    # https://www.researchgate.net/post/How_to_calculate_OR_odd_ratio_if_one_of_groups_is_0_in_a_case-control_study
    # Add 0.5 to each cell to avoid zero value cell problem
    correctionTerm<-0.5
    
    n11<-neighborList$localDegree[j] + correctionTerm
    n10<-neighborList$globalDegree[j] - neighborList$localDegree[j] + correctionTerm
    n01<-neighborList$numOfgeneOverlap - neighborList$localDegree[j] + correctionTerm
    n00<-neighborList$numOfgraphReducedGene - neighborList$globalDegree[j] - n01 + correctionTerm
    
    neighborList$oddsRatio[j]<-(log(n11) + log(n00) - log(n10) - log(n01))
  }
  
  neighborListFrame<-data.frame(neighborList$reducedGraphId,neighborList$name,neighborList$localDegree,neighborList$globalDegree,neighborList$pValueRaw,neighborList$oddsRatio,stringsAsFactors = FALSE)
  colnames(neighborListFrame)<-c("idx","name","localDegree","globalDegree","pValueRaw","oddsRatio")
  neighborListFrame<-neighborListFrame[order(neighborListFrame$pValueRaw),]
  
  localDegreeCutoff<-2
  neighborListFrame<-neighborListFrame[neighborListFrame$localDegree>=localDegreeCutoff,]
  
  cat(sprintf("Only test neighbor nodes with local degree equals or exceeds %s\n",localDegreeCutoff))
  cat(sprintf("Multiple hypothesis corrections for %s neighbor nodes in the network\n",nrow(neighborListFrame)))
  
  neighborListFrame$pValueFDR<-p.adjust(neighborListFrame$pValueRaw,method="BH")
  neighborListFrame$pValueBonferroni<-p.adjust(neighborListFrame$pValueRaw,method="bonferroni")
  
  #write.table(neighborListFrame,file="result.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
  #linkerListFrame<-neighborListFrame[neighborListFrame$pValueRaw < pValueCutoff, ] 
  
  if( pValueAdj == "BH" ) {
    linkerListFrame<-neighborListFrame[ neighborListFrame$pValueFDR < pValueCutoff, ] 
  }
  
  if( pValueAdj == "Bonferroni" ) {
    linkerListFrame<-neighborListFrame[ neighborListFrame$pValueBonferroni < pValueCutoff, ] 
  }
  
  cat(sprintf("For p-value %s cut-off, %s nodes were included as linker nodes\n",pValueCutoff,dim(linkerListFrame)[1]))
  linkerListFrame$type<-rep("linker",dim(linkerListFrame)[1])
  
  # include linker nodes into the list of initial candidate nodes 
  selectedGene<-{}
  selectedGene$idx<-union(geneOverlap$idx,linkerListFrame$idx)
  selectedGene$type<-c(geneOverlap$type,linkerListFrame$type)
  selectedGene$name<-c(geneOverlap$name,as.character(linkerListFrame$name))
  selectedNodeType<-data.frame(selectedGene$name,selectedGene$type,stringsAsFactors = FALSE)
  colnames(selectedNodeType)<-c("name","type") 
  
  
  
  cat(sprintf("Connecting %s candidate nodes and %s linker nodes\n",length(geneOverlap$name),dim(linkerListFrame)[1]))
  # sub-network of linker nodes and candidate nodes
  graphTemp<-induced.subgraph(graphReduced,selectedGene$idx)
  
  # add vertex attributes
  graphTemp<-set_vertex_attr(graphTemp,name="nodeType",index=V(graphTemp),value="candidate")
  graphTemp<-set_vertex_attr(graphTemp,name="nodeType",index=as.character(linkerListFrame$name),value="linker")
  
  
  # remove isolated vertex
  isolatedNodes<-names(which(degree(graphTemp)<1))
  
  if( keepIsolatedNodes ) {
    graphOutput<-graphTemp
    cat(sprintf("Keep %s isolated candidate nodes from the input\n",length(isolatedNodes)))
  } else {
    graphOutput<-delete.vertices(graphTemp,which(degree(graphTemp)<1))
    cat(sprintf("Remove %s isolated candidate nodes from the input\n",length(isolatedNodes)))
  }
  
  #numOfNodes<-length(V(graphTemp))
  #numOfEdges<-length(E(graphTemp))
  #cat(sprintf("Loading network of %s nodes and %s interactions\n",numOfNodes,numOfEdges))
  
  numOfNodes<-length(V(graphOutput))
  numOfEdges<-length(E(graphOutput))
  
  cat(sprintf("Final network contains %s nodes and %s interactions\n",numOfNodes,numOfEdges))
  
  communities <- list()
  # Assign community membership for the sub-network
  
  if (communityMethod == "ebc" )
  {
    cat(sprintf("Detecting modules using \"edge betweeness\" method\n"))
    community <- edge.betweenness.community(graphOutput)
    #communities$`Edge betweenness` <- ebc
    moduleMembership<-membership(community)
  }
  
  if (communityMethod == "lec" )
  {
    cat(sprintf("Detecting modules using \"leading eigenvector\" method\n"))
    community<- leading.eigenvector.community(graphOutput,options=list(maxiter=1000000))
    moduleMembership<-membership(community)
    
  }
  
  moduleMembership<-moduleMembership[names(moduleMembership)]
  moduleMembershipFrame<-data.frame(names(moduleMembership),moduleMembership,stringsAsFactors = FALSE)
  colnames(moduleMembershipFrame)<-c("geneSymbol","membership")
  moduleMembershipFrame<-moduleMembershipFrame[order(moduleMembershipFrame$membership),]
  
  
  #write.table(memb.ebc,file="memb.ebc.txt",sep="\t",quote=FALSE,col.names=FALSE)
  
  #graphOutput$layout<-layout.drl(graphOutput,options=list(simmer.attraction=0,edge.cut=1,expansion.attraction=0))
  #plot(ebc,graphOutput,vertex.size=3,vertex.label.dist=1.5)
  
  netbox<-{}
  #netbox$name<-V(graphOutput)$name
  # bug : get.edges only return half number of edges, use get.edgelist #
  #netbox$edgelist<-get.edgelist(graphOutput)
  #netbox$interactionType<-rep("INTERACT",length(netbox$edgelist[,1]))
  
  #netbox$edgelist<-get.edges(graphOutput,V(graphOutput))
  #netboxTemp<-data.frame(netbox$edgelist[,1],netbox$interactionType,netbox$edgelist[,2])
  #colnames(netboxTemp)<-c("geneAidx","interaction","geneBidx")
  
  #netboxTemp$geneSymbolA<-get.vertex.attribute(graphOutput,name="name",index=netboxTemp$geneAidx)  
  #netboxTemp$geneSymbolB<-get.vertex.attribute(graphOutput,name="name",index=netboxTemp$geneBidx)  
  
  #netboxOutput<-data.frame(netboxTemp$geneSymbolA,netboxTemp$interaction,netboxTemp$geneSymbolB)
  #colnames(netboxOutput)<-c("geneA","interaction","geneB")
  
  netbox$edgelist<-get.edgelist(graphOutput)
  
  ####
  edgeLabelList<-get.edge.attribute(graphOutput)$INTERACTION_TYPE
  mergedEdgeLabels<-lapply(edgeLabelList,function(x) 
  {   
    
    content<-unlist(unique(strsplit(x,";")))
    contentMerged<-paste(content,collapse=";")
    numOfcontent<-length(content)
    data.frame(numOfcontent,contentMerged,stringsAsFactors = FALSE)
    #paste(x,collapse=",")
  #},mc.cores=useCores)
  })
    
  mergedEdgeLabels<-do.call(rbind,mergedEdgeLabels)
  
  ####
  #netbox$interactionType<-rep("INTERACT",length(netbox$edgelist[,1]))
  
  #netboxOutput<-data.frame(netbox$edgelist[,1],netbox$interactionType,netbox$edgelist[,2],stringsAsFactors = FALSE)
  netboxOutput<-data.frame(netbox$edgelist[,1],mergedEdgeLabels$contentMerged,netbox$edgelist[,2],stringsAsFactors = FALSE)
  colnames(netboxOutput)<-c("geneA","interaction","geneB")
  
  if( keepIsolatedNodes ) {
    netboxOutputExtra<-data.frame(isolatedNodes,rep("INTERACT",length(isolatedNodes)),isolatedNodes,stringsAsFactors = FALSE)
    colnames(netboxOutputExtra)<-c("geneA","interaction","geneB")
    netboxOutput<-rbind(netboxOutput,netboxOutputExtra)
  }
  
  #write.table(netboxOutput,file="network.sif",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
  
  result<-list(netboxGraph=graphOutput,
               netboxCommunity=community,
               netboxOutput=netboxOutput,
               nodeType=selectedNodeType, 
               moduleMembership=moduleMembershipFrame, 
               neighborData=neighborListFrame)
  
  return(result)
}
