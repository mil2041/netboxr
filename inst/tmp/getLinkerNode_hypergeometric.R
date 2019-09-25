library(plyr)
library(netrankr)

iterations<-1000
useCores<-6

time_start<-proc.time()
neighborListFrameList<-list()   
originalNeighborListFrame<-neighborListFrame
  
  
neighborListFrameList<-mclapply(1:iterations, function(m)
{  
  
  graphReduced2<-rewire(graphReduced, with = keeping_degseq(niter = vcount(graphReduced) * 10))
  # Get the genes from the input that overlap with the vertexs in the network  
  geneOverlap<-{}  
  geneOverlap$name<-intersect(geneList,V(graphReduced2)$name)
  geneOverlap$idx<-match(geneOverlap$name,V(graphReduced2)$name)
  geneOverlap$type<-rep("candidate",length(geneOverlap$name))
  cat(sprintf("%s / %s candidate nodes match the name in the network of %s nodes \n",length(geneOverlap$name),length(geneList),length(V(graphReduced)$name)))
  
  geneOverlapFrame<-data.frame(geneOverlap$name,geneOverlap$idx,stringsAsFactors = FALSE)
  colnames(geneOverlapFrame)<-c("name","idx")
  
  #####
  #####
  
  neighborList<-{}
  neighborList$numOfgraphReducedGene<-length(V(graphReduced2))
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
    
    neighborTemp<-neighbors(graphReduced2,v=geneOverlap$idx[k],mode="all")
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
  
  neighborList$globalDegree<-degree(graphReduced2,v=neighborList$reducedGraphId,mode="all")
  
  
  ###
  
  out<-lapply( 1:length(neighborList$reducedGraphId), function(j)
  {  
    
    idx<-neighborList$reducedGraphId[j]
    name<-neighborList$name[j]
    globalNeighborId<-neighbors(graphReduced2,v=idx,mode=1)
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
  
  return(neighborListFrame)
  
},mc.cores=useCores)

cc<-rbind.fill(neighborListFrameList)

#####

pValueSample<-{}

for(selectedName in unique(originalNeighborListFrame$name)){
  
  observedScore<-(-1)*log10(originalNeighborListFrame[originalNeighborListFrame$name %in% selectedName,]$pValueRaw)
  simulatedScores<-(-1)*log10(cc[cc$name %in% selectedName,]$pValueRaw)
  pValueSample[selectedName]<-(sum(observedScore>=simulatedScores)+1) / (length(cc[cc$name %in% selectedName,]$pValueRaw) +1)
  #pValueSample[selectedName]<-(sum(originalNeighborListFrame[originalNeighborListFrame$name %in% selectedName,]$pValueRaw >= cc[cc$name %in% selectedName,]$pValueRaw)) / (length(cc[cc$name %in% selectedName,]$pValueRaw))
  
}

dd<-data.frame(names(pValueSample),pValueSample,stringsAsFactors = FALSE)

dd<-dd[order(dd$pValueSample),]
DT::datatable(dd, rownames = FALSE)

dd$pvalueFDR<-p.adjust(dd$pValueSample,method="BH")
DT::datatable(dd, rownames = FALSE)
nrow(dd[dd$pvalueFDR<=0.05,])
head(dd)
intersect(rownames(dd[dd$pvalueFDR<=0.05,]), linkerDF[linkerDF$pValueFDR<0.05,]$name)

ss<-dd[dd$pvalueFDR<=0.05,][,1]
#length(ss)

  time_elapesed<-proc.time()-time_start
  cat (sprintf ("Time for calculating p-value of %s candidates: %.2f sec\n", length(neighborListFrame$name), time_elapesed[3]) )
  
  
threshold<-0.05
ff<-data.frame(neighborListFrame$name,pRaw,stringsAsFactors = FALSE)
ff$fdr<-p.adjust(ff$pRaw,method="BH")
ff<-ff[order(ff$pRaw,ff$fdr),]
ff[ff$fdr<threshold,]
nrow(ff[ff$fdr<threshold,])
head(ff,20)

linkerNodeSelected<-ff[ff$fdr<threshold,]$neighborListFrame.name

