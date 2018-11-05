#' @title Simplify sif network into igraph network graph object with merged edge attributes 
#' 
#' @description  
#' This function removes duplicated edges and loops to create an igraph graph 
#' object from tab delimited sif formatted network file.  
#' 
#' @details 
#' For undirected graph, \code{networkSimplify} removes duplicated edges 
#' and loops to create an igraph graph object from tab delimited sif 
#' formatted network file.  
#' 
#' For directed graph, \code{networkSimplify} selects the first edge and
#' removes the rest duplicated edges and loops to create an igraph graph
#' object from tab delimited sif
#' formatted network file.       
#' 
#' @param sifNetwork A file with sif network format (There are three columns 
#'        in the file separated by tab, nodeA interactionType nodeB )  
#' 
#' @param directed Logical, treat network as directed or undirected graph
#' @param useCores number of cores 
#' 
#' @return a igraph graph object
#'   
#' @author Eric Minwei Liu, \email{emliu.research@gmail.com}
#' 
#' @examples
#' data(netbox2010)
#'
#' sifNetwork<-netbox2010$network         
#' graphReduced<-networkSimplify(sifNetwork,directed = FALSE)
#'
#' @concept netboxr
#' @export
#' @import igraph
#' @import parallel
networkSimplify2<-function(sifNetwork,directed=FALSE,useCores=1){
  
  time_start<-proc.time()
  
  relationsTable<-sifNetwork
  originalColumnName<-colnames(relationsTable)
  originalColumnName[1:6]<-c("PARTICIPANT_A","INTERACTION_TYPE","PARTICIPANT_B",
                             "INTERACTION_DATA_SOURCE","INTERACTION_PUBMED_ID","PATHWAY_NAMES")
  #relationsTable<-sifNetwork
  #relationsTable<-read.table("PCWithLoc_PPI_05092014_simple.sif",header=FALSE,sep="\t")
  removeColumn=c(1,3)
  relationsTableAttr<-relationsTable[,-removeColumn,drop=FALSE]
  
  
  interactions<-data.frame(relationsTable[, 1],relationsTable[, 3],relationsTableAttr,stringsAsFactors=FALSE)
  colnames(interactions)<-c(originalColumnName[1],originalColumnName[3],colnames(relationsTableAttr))
  
  # load graph as un-directed or directed graph
  graphFull <- graph.data.frame(interactions[,1:6], directed=directed)
  numOfNodes<-length(V(graphFull))
  numOfEdges<-length(E(graphFull))
  cat(sprintf("Loading network of %s nodes and %s interactions\n",numOfNodes,numOfEdges))
  if( directed == TRUE ){
    directionality<-"directed"  
  }else{
    directionality<-"undirected"  
  }  
  cat(sprintf("Treated as %s network \n",directionality))
  
  #print(graphFull, e=TRUE, v=TRUE)
  ## The opposite operation
  #get.data.frame(graphFull, what="vertices")
  #get.data.frame(graphFull, what="edges")
  
  # remove multiple interactions among the same pair of nodes and remove interaction loops
  graphReduced<-simplify(graphFull, remove.multiple=TRUE, remove.loops=TRUE, 
                         edge.attr.comb=c(INTERACTION_TYPE="concat", 
                                          INTERACTION_DATA_SOURCE="concat",
                                          INTERACTION_PUBMED_ID="concat",
                                          PATHWAY_NAMES="concat")
                         
  )
  
  #####
  #  merge edge attributes 1
  #####
  
  cat(sprintf("Merging edge attribute: INTERACTION_TYPE \n"))
  
  edgeLabelList<-get.edge.attribute(graphReduced)$INTERACTION_TYPE
  mergedEdgeLabels<-mclapply(edgeLabelList,function(x) 
  {   
    
    content<-unique(unlist(strsplit(x,";")))
    contentMerged<-paste(content,collapse=";")
    numOfcontent<-length(content)
    data.frame(numOfcontent,contentMerged,stringsAsFactors = FALSE)
    #paste(x,collapse=",")
  },mc.cores=useCores)
  
  mergedEdgeLabels<-do.call(rbind,mergedEdgeLabels)
  
  #mergedEdgeLabels<-unlist(mergedEdgeLabels)
  
  graphReduced<-set_edge_attr(graphReduced,"INTERACTION_TYPE",value=mergedEdgeLabels$contentMerged)
  graphReduced<-set_edge_attr(graphReduced,"INTERACTION_TYPE_COUNTS",value=mergedEdgeLabels$numOfcontent)
  
  
  #####
  #  merge edge attributes 2
  #####
  
  cat(sprintf("Merging edge attribute: INTERACTION_DATA_SOURCE \n"))
  
  edgeLabelList<-get.edge.attribute(graphReduced)$INTERACTION_DATA_SOURCE
  mergedEdgeLabels<-lapply(edgeLabelList,function(x) 
  {   
    
    content<-unique(unlist(strsplit(x,";")))
    contentMerged<-paste(content,collapse=";")
    numOfcontent<-length(content)
    data.frame(numOfcontent,contentMerged,stringsAsFactors = FALSE)
    #paste(x,collapse=",")
  })
  
  mergedEdgeLabels<-do.call(rbind,mergedEdgeLabels)
  
  #mergedEdgeLabels<-unlist(mergedEdgeLabels)
  
  graphReduced<-set_edge_attr(graphReduced,"INTERACTION_DATA_SOURCE",value=mergedEdgeLabels$contentMerged)
  graphReduced<-set_edge_attr(graphReduced,"INTERACTION_DATA_SOURCE_COUNTS",value=mergedEdgeLabels$numOfcontent)
  
  
  
  #####
  #  merge edge attributes 3
  #####
  
  cat(sprintf("Merging edge attribute: INTERACTION_PUBMED_ID \n"))
  
  edgeLabelList<-get.edge.attribute(graphReduced)$INTERACTION_PUBMED_ID
  mergedEdgeLabels<-lapply(edgeLabelList,function(x) 
  {   
    
    content<-unique(unlist(strsplit(x,";")))
    contentMerged<-paste(content,collapse=";")
    numOfcontent<-length(content)
    data.frame(numOfcontent,contentMerged,stringsAsFactors = FALSE)
    #paste(x,collapse=",")
  })
  
  mergedEdgeLabels<-do.call(rbind,mergedEdgeLabels)
  
  #mergedEdgeLabels<-unlist(mergedEdgeLabels)
  
  graphReduced<-set_edge_attr(graphReduced,"INTERACTION_PUBMED_ID",value=mergedEdgeLabels$contentMerged)
  graphReduced<-set_edge_attr(graphReduced,"INTERACTION_PUBMED_ID_COUNTS",value=mergedEdgeLabels$numOfcontent)
  
  #####
  #  merge edge attributes 4
  #####
  
  cat(sprintf("Merging edge attribute: PATHWAY_NAMES \n"))
  
  edgeLabelList<-get.edge.attribute(graphReduced)$PATHWAY_NAMES
  mergedEdgeLabels<-lapply(edgeLabelList,function(x) 
  {   
    
    content<-unique(unlist(strsplit(x,";")))
    contentMerged<-paste(content,collapse=";")
    numOfcontent<-length(content)
    data.frame(numOfcontent,contentMerged,stringsAsFactors = FALSE)
    #paste(x,collapse=",")
  })
  
  mergedEdgeLabels<-do.call(rbind,mergedEdgeLabels)
  
  #mergedEdgeLabels<-unlist(mergedEdgeLabels)
  
  graphReduced<-set_edge_attr(graphReduced,"PATHWAY_NAMES",value=mergedEdgeLabels$contentMerged)
  graphReduced<-set_edge_attr(graphReduced,"PATHWAY_NAMES_COUNTS",value=mergedEdgeLabels$numOfcontent)
  
  #####
  
  numOfNodes<-length(V(graphReduced))
  numOfEdges<-length(E(graphReduced))
  cat(sprintf("Removing multiple interactions and loops\n"))
  cat(sprintf("Returning network of %s nodes and %s interactions\n",numOfNodes,numOfEdges))
  
  time_elapesed<-proc.time()-time_start
  cat (sprintf ("Time for network simplify: %.2f sec\n", time_elapesed[3]) )
  
  
  return(graphReduced)
   
}
