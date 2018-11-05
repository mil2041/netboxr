library(plyr)
library(netrankr)

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

#####

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

