library(plyr)
library(netrankr)

iterations<-1000
useCores<-6

time_start<-proc.time()

vertex_betweenness_list<-mclapply(1:iterations, function(m)
{  
  graphReduced2<-rewire(graphReduced, with = keeping_degseq(niter = vcount(graphReduced) * 10))
  #graphTemp2<-induced.subgraph(graphReduced2,union(neighborList$name[i], geneOverlap$name))
  #vertex_betweenness_vector2<-betweenness(graphTemp2,v=neighborList$name[i],directed=FALSE,normalized = TRUE)
  
  graphTemp2<-induced.subgraph(graphReduced2,union(neighborListFrame$name, geneOverlap$name))
  #graphTemp3<-delete.vertices(graphTemp2,which(degree(graphTemp2)<1))
  #vertex_betweenness_vector2<-betweenness(graphTemp2,vids=neighborListFrame$name,mode="all",normalized = TRUE)
  #vertex_degree_vector2<-degree(graphTemp2,v=neighborListFrame$name,mode="all",normalized = TRUE)
  vertex_betweenness_vector2<-betweenness(graphTemp2,v=neighborListFrame$name,directed=FALSE,normalized = TRUE)
  
  #vertex_betweenness_vector2<-graphTemp3 %>% 
  #                        indirect_relations(type="depend_sp") %>% 
  #                        aggregate_positions(type="sum")
  
  return(vertex_betweenness_vector2)
  
},mc.cores=useCores)

s1<-unlist(vertex_betweenness_list)

#####

pRaw<-{}


  graphTemp<-induced.subgraph(graphReduced,union(neighborListFrame$name, geneOverlap$name))
  #graphTemp2<-delete.vertices(graphTemp,which(degree(graphTemp)<1))
  vertex_betweenness_vector<-betweenness(graphTemp,v=neighborListFrame$name,directed=FALSE,normalized = TRUE)
  #vertex_betweenness_vector2<-estimate_betweenness(graphTemp2,v=neighborListFrame$name,directed=FALSE,cutoff=3)
  #vertex_degree_vector<-degree(graphTemp,v=neighborListFrame$name[i],mode="all",normalized = TRUE)
  
  #start_time <- Sys.time()
  #vertex_betweenness_vector2<-estimate_betweenness(graphTemp2,v=neighborListFrame$name,directed=FALSE,cutoff=3)
  #end_time<-Sys.time()
  #end_time - start_time
  
  #start_time <- Sys.time()
  #vertex_betweenness_vector0<-betweenness(graphTemp2,v=neighborListFrame$name,directed=FALSE,normalized = FALSE)
  #end_time<-Sys.time()
  #end_time - start_time
  
  #start_time <- Sys.time()
  #vertex_betweenness_vector<-graphTemp2 %>% 
  #  indirect_relations(type="depend_sp") %>% 
  #  aggregate_positions(type="sum")

  #end_time<-Sys.time()
  #end_time - start_time
  
for(i in 1:length(neighborListFrame$name)){
  
  vertex_betweenness_vector0<-vertex_betweenness_vector[names(vertex_betweenness_vector)==neighborListFrame$name[i]]
  
  s2<-s1[names(s1) %in% neighborListFrame$name[i]]
  
  pRaw[i]<-(sum(s2>=vertex_betweenness_vector0)+1)/(iterations+1)
  #pRaw[i]<-(sum(s2>=vertex_degree_vector)+1)/(iterations+1)
  
  cat(sprintf("%s / %s name:%s : pvalue: %s \n",i,length(neighborListFrame$name),neighborListFrame$name[i],pRaw[i] ))
  
  
}

  time_elapesed<-proc.time()-time_start
  cat (sprintf ("Time for calculating p-value of %s candidates: %.2f sec\n", length(neighborListFrame$name), time_elapesed[3]) )
  
  
threshold<-0.05
ff<-data.frame(neighborListFrame$name,pRaw,stringsAsFactors = FALSE)
ff$fdr<-p.adjust(ff$pRaw,method="BH")
ff<-ff[order(ff$pRaw,ff$fdr),]
head(ff[ff$fdr<threshold,])
nrow(ff[ff$fdr<threshold,])
head(ff,20)

#colnames(ff)<-c("name","pRaw","pFDR")
linkerNodeSelected<-ff[ff$fdr<threshold,]$neighborListFrame.name



