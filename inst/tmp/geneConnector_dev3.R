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
#' @import plyr rbind.fill
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
  neighborList$oddsRatio<-{}
  
  #for (i in 1:length(neighborList$reducedGraphId))
  #{
  #  neighborList$name[i]<-get.vertex.attribute(graphReduced,name="name",index=neighborList$reducedGraphId[i])  
  #}  
  
  tmpFrame<-neighborTotalFrame[neighborTotalFrame$idx %in% neighborList$reducedGraphId,]
  rownames(tmpFrame)<-tmpFrame$idx
  tmpFrame<-tmpFrame[match(neighborList$reducedGraphId,rownames(tmpFrame)),]
  
  neighborList$name<-tmpFrame$name
  
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
    
    localDegree<-a
    globalDegree<-b
    pValueRaw<-neighborList$pValueRaw[j]
    oddsRatio<-neighborList$oddsRatio[j]
    
    out<-data.frame(idx,name,localDegree,globalDegree,pValueRaw,oddsRatio,stringsAsFactors = FALSE)
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
  
  #######
  
  
  
}
