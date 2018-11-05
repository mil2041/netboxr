#' Generate sub-network mapping from a candidate gene list version 2
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
geneConnector2<-function(geneList,networkGraph,directed=FALSE,pValueAdj="BH",pValueCutoff=0.05,communityMethod="lec",keepIsolatedNodes=FALSE){
  
  graphReduced<-networkGraph
  
  # Get the genes from the input that overlap with the vertexs in the network  
  geneOverlap<-{}  
  geneOverlap$name<-intersect(geneList,V(graphReduced)$name)
  geneOverlap$idx<-match(geneOverlap$name,V(graphReduced)$name)
  geneOverlap$type<-rep("candidate",length(geneOverlap$name))
  cat(sprintf("%s / %s candidate nodes match the name in the network of %s nodes \n",length(geneOverlap$name),length(geneList),length(V(graphReduced)$name)))
  
  geneOverlapFrame<-data.frame(geneOverlap$name,geneOverlap$idx,stringsAsFactors = FALSE)
  colnames(geneOverlapFrame)<-c("name","idx")
  
  #####
  #####
  
  neighborList<-{}
  neighborList$numOfgraphReducedGene<-length(V(graphReduced))
  neighborList$numOfgeneOverlap<-length(geneOverlap$idx)  
  
  #neighborList$reducedGraphId<-neighbors(graphReduced,v=geneOverlap$idx,mode=1)
  
  neighborTotal<-{}
  neighborTemp<-{}
  
  #for(k in 1:length(geneOverlap$idx))
  #{
  #  neighborTemp<-neighbors(graphReduced,v=geneOverlap$idx[k],mode="all")
  #  neighborTotal<-union(neighborTotal,neighborTemp)
  #}
  
  neighborTotal<-lapply(1:length(geneOverlap$idx), function(k) 
  {
  
    neighborTemp<-neighbors(graphReduced,v=geneOverlap$idx[k],mode="all")
    #neighborTotal<-union(neighborTotal,neighborTemp)
    return(neighborTemp)
    
  })
  
  neighborsVector<-unlist(neighborTotal)
  neighborTotalFrame<-data.frame(names(neighborsVector),neighborsVector,stringsAsFactors = FALSE)
  colnames(neighborTotalFrame)<-c("name","idx")
  
  # This frame may contains the original gene nodes, 
  neighborTotalFrame<-neighborTotalFrame[(!duplicated(neighborTotalFrame$idx)),]
  
  neighborList$reducedGraphId<-setdiff(neighborTotalFrame$idx, geneOverlap$idx)
  
  
  ######
  
  neighborList$name<-{}
  neighborList$localDegree<-{}
  
  #for (i in 1:length(neighborList$reducedGraphId))
  #{
  #  neighborList$name[i]<-get.vertex.attribute(graphReduced,name="name",index=neighborList$reducedGraphId[i])  
  #}  
  
  tmpFrame<-neighborTotalFrame[neighborTotalFrame$idx %in% neighborList$reducedGraphId,]
  rownames(tmpFrame)<-tmpFrame$idx
  tmpFrame<-tmpFrame[match(neighborList$reducedGraphId,rownames(tmpFrame)),]
  
  neighborList$name<-tmpFrame$name
  
  
  #####
  
  graphTemp<-induced.subgraph(graphReduced,union(neighborList$reducedGraphId[i], geneOverlap$idx))
  
  vertex_betweenness_vector<-betweenness(graphTemp,v=neighborList$name[i],directed=FALSE,normalized = TRUE)
  
  vertex_betweenness_frame<-data.frame(names(vertex_betweenness_vector),vertex_betweenness_vector,stringsAsFactors = FALSE)
  colnames(vertex_betweenness_frame)<-c("name","vertexBetweeness")
  vertex_betweenness_frame<-vertex_betweenness_frame[order(-vertex_betweenness_frame$vertexBetweeness),]
  vertex_betweenness_frame2<-vertex_betweenness_frame[vertex_betweenness_frame$name %in% neighborList$name,]
  
  #####
  
  graphTemp<-induced.subgraph(graphReduced,union(neighborList$name, geneOverlap$name))
  
  vertex_betweenness_vector<-betweenness(graphTemp,v=neighborList$name,directed=FALSE,normalized = TRUE)
  
  vertex_betweenness_frame<-data.frame(names(vertex_betweenness_vector),vertex_betweenness_vector,stringsAsFactors = FALSE)
  colnames(vertex_betweenness_frame)<-c("name","vertexBetweeness")
  vertex_betweenness_frame<-vertex_betweenness_frame[order(-vertex_betweenness_frame$vertexBetweeness),]
  vertex_betweenness_frame2<-vertex_betweenness_frame[vertex_betweenness_frame$name %in% neighborList$name,]
  
  ######
  
  iterations<-1000
  useCores<-6
  
  pRaw<-{}
  
  for(i in 1:length(neighborListFrame$name)){
  
  graphTemp<-induced.subgraph(graphReduced,union(neighborListFrame$name, geneOverlap$name))
  
  vertex_betweenness_vector<-betweenness(graphTemp,v=neighborListFrame$name[i],directed=FALSE,normalized = TRUE)
  
    
    
    vertex_betweenness_list<-mclapply(1:iterations, function(m)
    {  
      graphReduced2<-rewire(graphReduced, with = keeping_degseq(niter = vcount(graphReduced) * 10))
      #graphTemp2<-induced.subgraph(graphReduced2,union(neighborList$name[i], geneOverlap$name))
      #vertex_betweenness_vector2<-betweenness(graphTemp2,v=neighborList$name[i],directed=FALSE,normalized = TRUE)
    
      graphTemp2<-induced.subgraph(graphReduced2,union(neighborListFrame$name, geneOverlap$name))
      vertex_betweenness_vector2<-betweenness(graphTemp2,v=neighborListFrame$name,directed=FALSE,normalized = TRUE)
      return(vertex_betweenness_vector2)
      
    },mc.cores=useCores)
    
    s1<-unlist(vertex_betweenness_list)
    
    
    
    
    #sum(s1>=vertex_betweenness_vector)/iterations
    
    s2<-s1[names(s1) %in% neighborListFrame$name[i]]
    
    pRaw[i]<-sum(s2>=vertex_betweenness_vector)/iterations
  
    cat(sprintf("%s / %s name:%s : pvalue: %s \n",i,length(neighborListFrame$name),neighborListFrame$name[i],pRaw[i] ))
    
  }  
    
}  

#######
iterations<-5000
useCores<-6

vertex_degree_list<-mclapply(1:iterations, function(m)
{  
  graphReduced2<-rewire(graphReduced, with = keeping_degseq(niter = vcount(graphReduced) * 10))
  #graphTemp2<-induced.subgraph(graphReduced2,union(neighborList$name[i], geneOverlap$name))
  #vertex_betweenness_vector2<-betweenness(graphTemp2,v=neighborList$name[i],directed=FALSE,normalized = TRUE)
  
  graphTemp2<-induced.subgraph(graphReduced2,union(neighborListFrame$name, geneOverlap$name))
  #vertex_betweenness_vector2<-closeness(graphTemp2,vids=neighborListFrame$name,mode="all",normalized = TRUE)
  vertex_degree_vector2<-degree(graphTemp2,v=neighborListFrame$name,mode="all",normalized = TRUE)
  
  return(vertex_degree_vector2)
  
},mc.cores=useCores)

s1<-unlist(vertex_degree_list)

#######

i<-327

pRaw<-{}

for(i in 1:length(neighborListFrame$name)){
  
  graphTemp<-induced.subgraph(graphReduced,union(neighborListFrame$name, geneOverlap$name))
  #graphTemp2<-delete.vertices(graphTemp,which(degree(graphTemp)<1))
  vertex_betweenness_vector<-betweenness(graphTemp,v=neighborListFrame$name[i],directed=FALSE,normalized = TRUE)
  #vertex_degree_vector<-degree(graphTemp,v=neighborListFrame$name[i],mode="all",normalized = TRUE)
  
  s2<-s1[names(s1) %in% neighborListFrame$name[i]]
  
  pRaw[i]<-(sum(s2>=vertex_betweenness_vector)+1)/(iterations+1)
  #pRaw[i]<-(sum(s2>=vertex_degree_vector)+1)/(iterations+1)
  
  cat(sprintf("%s / %s name:%s : pvalue: %s \n",i,length(neighborListFrame$name),neighborListFrame$name[i],pRaw[i] ))
  
    
}

 threshold<-0.05
 ff<-data.frame(neighborListFrame$name,pRaw,stringsAsFactors = FALSE)
 ff$fdr<-p.adjust(ff$pRaw,method="BH")
 ff<-ff[order(ff$fdr),]
 ff[ff$fdr<threshold,]
 
 linkerNodeSelected<-ff[ff$fdr<threshold,]$neighborListFrame.name
 
  ######
  
  neighborList$globalDegree<-degree(graphReduced,v=neighborList$reducedGraphId,mode="all")
  
  
  ###
  
  out<-lapply( 1:length(neighborList$reducedGraphId), function(j)
  {  
    
    idx<-neighborList$reducedGraphId[j]
    name<-neighborList$name[j]
    globalNeighborId<-neighbors(graphReduced,v=idx,mode=1)
    localNeighborId<-intersect(geneOverlap$idx,globalNeighborId)
    neighborList$localDegree[j]<-length(localNeighborId)
    #neighborList_localDegree<-length(localNeighborId)
    
    # For each linker node candidate (neighbor node)
    # a: (local degree) number of connections to nodes in reduced graph that overlaps with gene list.   
    # b: (global degree) number of connections to all the nodes in the reduced graph. 
    # c: number of "not" connected nodes (total number of nodes in the reduced graph minus number of conncted nodes)
    # d: number of nodes in the reduced graph that overlaps with gene list. 
    
    a<-neighborList$localDegree[j]
    b<-neighborList$globalDegree[j]
    c<-(neighborList$numOfgraphReducedGene - neighborList$globalDegree[j])  
    d<-neighborList$numOfgeneOverlap
    
    # hypergeometric distribution for p-value calculation
    #neighborList$pValueRaw[j]<-(1-phyper((a-1),b,c,d))
    
    # http://pedagogix-tagc.univ-mrs.fr/courses/ASG1/practicals/go_statistics_td/go_statistics_td_2015.html
    neighborList$pValueRaw[j]<-phyper((a-1),b,c,d, lower.tail = FALSE)
    
    localDegree<-a
    globalDegree<-b
    pValueRaw<-neighborList$pValueRaw[j]
    
    out<-data.frame(idx,name,localDegree,globalDegree,pValueRaw,stringsAsFactors = FALSE)
    return(out)
    
  })
  
  neighborListFrame<-rbind.fill(out)
  
  #neighborListFrame<-data.frame(neighborList$reducedGraphId,neighborList$name,neighborList$localDegree,neighborList$globalDegree,neighborList$pValueRaw,stringsAsFactors = FALSE)
  #colnames(neighborListFrame)<-c("idx","name","localDegree","globalDegree","pValueRaw")
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
  
  #
  linkerListFrame<-neighborListFrame[neighborListFrame$name %in% linkerNodeSelected,]
  
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
  netbox$interactionType<-rep("INTERACT",length(netbox$edgelist[,1]))
  
  netboxOutput<-data.frame(netbox$edgelist[,1],netbox$interactionType,netbox$edgelist[,2],stringsAsFactors = FALSE)
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
