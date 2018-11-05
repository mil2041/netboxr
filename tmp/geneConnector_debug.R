


#####

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
      
      a<-neighborList$localDegree[j]
      b<-neighborList$globalDegree[j]
      c<-(neighborList$numOfgraphReducedGene - neighborList$globalDegree[j])  
      d<-neighborList$numOfgeneOverlap
      
      # hypergeometric distribution for p-value calculation
      neighborList$pValueRaw[j]<-(1-phyper((a-1),b,c,d))
    }
    
    neighborListFrame<-data.frame(neighborList$reducedGraphId,neighborList$name,neighborList$localDegree,neighborList$globalDegree,neighborList$pValueRaw,stringsAsFactors = FALSE)
    colnames(neighborListFrame)<-c("idx","name","localDegree","globalDegree","pValueRaw")
    neighborListFrame<-neighborListFrame[order(neighborListFrame$pValueRaw),]
    
    localDegreeCutoff<-2
    neighborListFrame<-neighborListFrame[neighborListFrame$localDegree>=localDegreeCutoff,]
    
    cat(sprintf("Only test neighbor nodes with local degree equals or exceeds %s\n",localDegreeCutoff))
    cat(sprintf("Multiple hypothesis corrections for %s neighbor nodes in the network\n",nrow(neighborListFrame)))
    
    neighborListFrame$pValueFDR<-p.adjust(neighborListFrame$pValueRaw,method="BH")
    neighborListFrame$pValueBonferroni<-p.adjust(neighborListFrame$pValueRaw,method="bonferroni")


    originalNeighborListFrame<-neighborListFrame

#####

iterations<-1000

neighborListFrameList<-list()    
useCores<-6

    
#for(iter in 1:iterations){
neighborListFrameList<-mclapply( 1:iterations, function(iter) {

    netboxGraph<-graphReduced
    
    #graphReduced<-networkGraph
    
    graphTmp<-rewire(netboxGraph, with = keeping_degseq(niter = vcount(netboxGraph) * 10))
  
    # Get the genes from the input that overlap with the vertexs in the network  
    geneOverlap<-{}  
    geneOverlap$name<-intersect(geneList,V(graphTmp)$name)
    geneOverlap$idx<-match(geneOverlap$name,V(graphTmp)$name)
    geneOverlap$type<-rep("candidate",length(geneOverlap$name))
    cat(sprintf("%s / %s candidate nodes match the name in the network of %s nodes \n",length(geneOverlap$name),length(geneList),length(V(graphTmp)$name)))
    
    neighborList<-{}
    neighborList$numOfgraphTmpGene<-length(V(graphTmp))
    neighborList$numOfgeneOverlap<-length(geneOverlap$idx)  
    
    #neighborList$reducedGraphId<-neighbors(graphTmp,v=geneOverlap$idx,mode=1)
    
    neighborTotal<-{}
    neighborTemp<-{}
    for (k in 1:length(geneOverlap$idx))
    {
      neighborTemp<-neighbors(graphTmp,v=geneOverlap$idx[k],mode="all")
      neighborTotal<-union(neighborTotal,neighborTemp)
    }
    neighborList$reducedGraphId<-setdiff(neighborTotal, geneOverlap$idx)
    
    
    neighborList$name<-{}
    neighborList$localDegree<-{}
    
    for (i in 1:length(neighborList$reducedGraphId))
    {
      neighborList$name[i]<-get.vertex.attribute(graphTmp,name="name",index=neighborList$reducedGraphId[i])  
    }  
    
    neighborList$globalDegree<-degree(graphTmp,v=neighborList$reducedGraphId,mode="all")
    
    for (j in 1:length(neighborList$reducedGraphId))
    {  
      
      globalNeighborId<-neighbors(graphTmp,v=neighborList$reducedGraphId[j],mode=1)
      localNeighborId<-intersect(geneOverlap$id,globalNeighborId)
      neighborList$localDegree[j]<-length(localNeighborId)
      
      a<-neighborList$localDegree[j]
      b<-neighborList$globalDegree[j]
      c<-(neighborList$numOfgraphTmpGene - neighborList$globalDegree[j])  
      d<-neighborList$numOfgeneOverlap
      
      # hypergeometric distribution for p-value calculation
      neighborList$pValueRaw[j]<-(1-phyper((a-1),b,c,d))
    }
    
    neighborListFrame<-data.frame(neighborList$reducedGraphId,neighborList$name,neighborList$localDegree,neighborList$globalDegree,neighborList$pValueRaw,stringsAsFactors = FALSE)
    colnames(neighborListFrame)<-c("idx","name","localDegree","globalDegree","pValueRaw")
    neighborListFrame<-neighborListFrame[order(neighborListFrame$pValueRaw),]
    
    localDegreeCutoff<-2
    neighborListFrame<-neighborListFrame[neighborListFrame$localDegree>=localDegreeCutoff,]
    
    cat(sprintf("Only test neighbor nodes with local degree equals or exceeds %s\n",localDegreeCutoff))
    cat(sprintf("Multiple hypothesis corrections for %s neighbor nodes in the network\n",nrow(neighborListFrame)))
    
    #neighborListFrame$pValueFDR<-p.adjust(neighborListFrame$pValueRaw,method="BH")
    #neighborListFrame$pValueBonferroni<-p.adjust(neighborListFrame$pValueRaw,method="bonferroni")

    #neighborListFrameList[[iter]]<-neighborListFrame
    
    return(neighborListFrame)
},mc.cores=useCores)    
    
cc<-rbind.fill(neighborListFrameList)

hist(-log10(cc[cc$name %in% "CDK6",]$pValueRaw))

originalNeighborListFrame[originalNeighborListFrame$name %in% "CDK6",]$pValueRaw

pValueSample<-{}

for(selectedName in unique(originalNeighborListFrame$name)){

pValueSample[selectedName]<-(sum(originalNeighborListFrame[originalNeighborListFrame$name %in% selectedName,]$pValueRaw >= cc[cc$name %in% selectedName,]$pValueRaw)+1) / (length(cc[cc$name %in% selectedName,]$pValueRaw) +1)
#pValueSample[selectedName]<-(sum(originalNeighborListFrame[originalNeighborListFrame$name %in% selectedName,]$pValueRaw >= cc[cc$name %in% selectedName,]$pValueRaw)) / (length(cc[cc$name %in% selectedName,]$pValueRaw))

}

dd<-data.frame(names(pValueSample),pValueSample,stringsAsFactors = FALSE)

dd<-dd[order(dd$pValueSample),]
DT::datatable(dd, rownames = FALSE)

dd$pvalueFDR<-p.adjust(dd$pValueSample,method="BH")
DT::datatable(dd, rownames = FALSE)
nrow(dd[dd$pvalueFDR<=0.05,])

intersect(rownames(dd[dd$pvalueFDR<=0.05,]), linkerDF[linkerDF$pValueFDR<0.05,]$name)

ss<-dd[dd$pvalueFDR<=0.05,][,1]

linkerListFrame<-originalNeighborListFrame[ originalNeighborListFrame$name %in% ss,]

#####

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



