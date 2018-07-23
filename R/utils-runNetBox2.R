#' Run NetBox2 network analysis
#' 
#' @param dataDir Where to store data 
#' @param cancerType study name from TCGA
#' @param mutationList User provided mutated gene list
#' @param ampGeneList User provided copy number amp alteration gene list
#' @param delGeneList User provided copy number deleted alteration gene list
#' @param epiSilencedList User provided epi genetic silencing gene list
#' @param mutationFreq  User provided mutated gene freq
#' @param ampGeneFreq User provided copy number amp alteration gene freq
#' @param delGeneFreq User provided copy number deleted alteration gene freq
#' @param epiSilencedFreq User provided epi genetic silencing gene freq
#' @param pathwayCommonsDb c("signal", "BioGrid")
#' @param directed TRUE or FALSE
#' @param linkerPValThreshold Default is 0.05
#' @param communityDetectionMethod "ebc", "lec", "none"
#' @param keepIsolatedNodes TRUE or FALSE
#' @param verbose Default is TRUE
#' 
#' @return null
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' selectedSampleId<-NA
#' #worDir<-getwd()
#' mutSig2CVthreshold<-0.1
#' rareMutationUpperLimit<-0.3
#' rareMutationLowerLimit<-0.1
#' rareMutationFreq<-0.02
#' 
#' #runNetBox2(dataDir,cancerType,
#' #           mutationList,ampGeneList,delGeneList,epiSilencedList,
#' #           mutationFreq,ampGeneFreq,delGeneFreq,epiSilencedFreq,
#' #           pathwayCommonsDb,directed,
#' #           linkerPValThreshold,communityDetectionMethod,
#' #           keepIsolatedNodes,verbose=TRUE)
#'
#' @concept netboxr
#' @export 
#' @importFrom HGNChelper checkGeneSymbols
#' @import gplots
#' @importFrom R.utils gzip
#' @importFrom plyr rbind.fill
runNetBox2<-function(dataDir,cancerType,
                     mutationList,ampGeneList,delGeneList,epiSilencedList,
                     mutationFreq,ampGeneFreq,delGeneFreq,epiSilencedFreq,
                     pathwayCommonsDb,directed,
                     linkerPValThreshold,communityDetectionMethod,
                     keepIsolatedNodes,verbose=TRUE){

  ampGeneListFiltered<-ampGeneList
  delGeneListFiltered<-delGeneList

  venn(list(ampCN=ampGeneListFiltered,delCN=delGeneListFiltered))

  fileName<-paste(cancerType,"_venn_geneList.png",sep="")
  filePath<-file.path(dataDir,cancerType)
  fileName<-file.path(filePath,fileName)
  png(file=fileName,width=1024,height=780)
  #venn(list(mutation=sigGenesMutation$gene,ampCN=ampGeneListFiltered,delCN=delGeneListFiltered))
  venn(list(mutation=mutationList,ampCN=ampGeneListFiltered,delCN=delGeneListFiltered,epiSilenced=epiSilencedList))
  dev.off()

  #intersect(epiSilencedList,ampGeneListFiltered)
  #intersect(mutationList,delGeneListFiltered)

  copyNumberGeneListFiltered<-union(ampGeneListFiltered,delGeneListFiltered)

  checkTable<-suppressWarnings(checkGeneSymbols(as.character(mutationList)))
  checkTable$freq<-mutationFreq

# convert to latest HGNC geneSymbol and remove obselete geneSymbol
  checkTable<-checkTable[!is.na(checkTable[,3]),]
  mutationList<-checkTable[,3]
  mutationFreq<-checkTable[,4]
  mutationDF<-data.frame(mutationList,mutationFreq,rep("Mutation",length(mutationList)))
  colnames(mutationDF)<-c("geneSymbol","freq","type")

####  
  
  if(length(ampGeneListFiltered)>0){
    checkTable<-suppressWarnings(checkGeneSymbols(as.character(ampGeneListFiltered)))
    checkTable$freq<-ampGeneFreq
# convert to latest HGNC geneSymbol and remove obselete geneSymbol
    checkTable<-checkTable[!is.na(checkTable[,3]),]
    ampGeneListFiltered<-checkTable[,3]
    ampGeneFreq<-checkTable[,4]
    ampDF<-data.frame(ampGeneListFiltered,ampGeneFreq,rep("Amp",length(ampGeneListFiltered)))
    colnames(ampDF)<-c("geneSymbol","freq","type")
  }else{ ampDF<-{}}

  if(length(delGeneListFiltered)>0){
    checkTable<-suppressWarnings(checkGeneSymbols(as.character(delGeneListFiltered)))
    checkTable$freq<-delGeneFreq

    # convert to latest HGNC geneSymbol and remove obselete geneSymbol
    checkTable<-checkTable[!is.na(checkTable[,3]),]
    delGeneListFiltered<-checkTable[,3]
    delGeneFreq<-checkTable[,4]
    delDF<-data.frame(delGeneListFiltered,delGeneFreq,rep("homDel",length(delGeneListFiltered)))
    colnames(delDF)<-c("geneSymbol","freq","type")
  }else{ delDF<-{}}


  alterationFreq<-rbind(mutationDF,ampDF,delDF)
  alterationFreq<-alterationFreq[order(alterationFreq$geneSymbol,-alterationFreq$freq),]
  #keep<-!duplicated(alterationFreq$geneSymbol)
  #alterationFreq<-alterationFreq[keep,]

  if(length(epiSilencedList)>0){
    hyperMetDF<-data.frame(epiSilencedList,rep(40,length(epiSilencedList)),rep("hyperMet",length(epiSilencedList)))
    colnames(hyperMetDF)<-c("geneSymbol","freq","type")
  }else{
    hyperMetDF<-{}
  }
    
#####

  ## check 573 cancer census genes
  #filePath<-file.path(PROJHOME,"annotation/Cosmic_census_list")
  #fileName<-"cosmic_CGC_0128_2016.txt"
  #fileName<-file.path(filePath,fileName)
  
  fileName<-system.file("extdata","cosmic_CGC_0128_2016.txt",package="netboxr")
  tmpData<-read.table(fileName,sep="\t",header=TRUE,
                    stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')

  ## check 720 epi enzymes
  #filePath<-file.path(PROJHOME,"annotation/epigenetic_gene_list/epifactors")
  #fileName<-"EpiGenes_main_genes.txt"
  #fileName<-file.path(filePath,fileName)
  
  fileName<-system.file("extdata","EpiGenes_main_genes.txt",package="netboxr")
  tmpData2<-read.table(fileName,sep="\t",header=TRUE,
                     stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')

  if( !is.null(hyperMetDF) ){
    alterationList<-rbind(mutationDF,ampDF,delDF,hyperMetDF)
  }else{
    alterationList<-rbind(mutationDF,ampDF,delDF)
  }
    
  alterationList$CGC<-rep("NA",nrow(alterationList))
  alterationList$epiEnzyme<-rep("NA",nrow(alterationList))

  alterationList[alterationList$geneSymbol %in% tmpData$Gene.Symbol,]$CGC<-"YES"
  alterationList[alterationList$geneSymbol %in% tmpData2$HGNC_symbol,]$epiEnzyme<-"YES"
  
  #####  

  #fusionList<-c("ERG","ETV1","ETV4","FLI1")

  #####
  
  copyNumberGeneListFiltered<-union(ampGeneListFiltered,delGeneListFiltered)
  copyNumberList<-copyNumberGeneListFiltered
  geneList<-union(mutationList,copyNumberList)
  geneList<-union(geneList,epiSilencedList)
  #geneList<-union(geneList,fusionList)
  length(geneList)

  ###
  #load network from internal pathwaycommons network
  
  networkDF<-{}
  
  if( "signal" %in% pathwayCommonsDb ){
  
    data(PC2V7_02212016_sif)
    network<-PC2V7_02212016_sif
    #dim(network)
    #names(network)
    networkDF[["signal"]]<-network
  }
  
  if( "BioGrid" %in% pathwayCommonsDb ){
    data(bioGridNetworkHomoSapiens34133)
    networkBioGrid<-bioGridNetworkHomoSapiens34133
    #dim(networkBioGrid)
    #names(networkBioGrid)

    networkBioGridE2<-networkBioGrid[networkBioGrid$EVIDENCE_NUM>=2,]
    networkDF[["BioGrid"]]<-networkBioGridE2[,1:5]
    
    #dim(networkBioGridE2)

  }
  
  selectedNetwork<-rbind.fill(networkDF)
  
  #networkBioGridE3<-networkBioGrid[networkBioGrid$EVIDENCE_NUM==1,]
  #trend<-table(networkBioGridE3$YEAR)
  #plot(trend)
  #networkBioGridE5<-networkBioGridE3[networkBioGridE3$YEAR>=2010,]
  #dim(networkBioGridE5)

  #networkBioGridE6<-rbind(networkBioGridE2,networkBioGridE5)

  #fileName<-"netbox2Network_w_Biogrid_e2.sif"
  fileName<-"selectedNetwork.sif.gz"
  filePath<-file.path(dataDir)
  fileName<-file.path(filePath,fileName)
  gz1<-gzfile(fileName,"w")
  write.table(selectedNetwork,file=gz1,sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
  close(gz1)
  
  #selectedNetwork2<-rbind(network,networkBioGridE6[,1:5])

##
  sifNetwork<-selectedNetwork
  graphReduced<-networkSimplify(sifNetwork,directed = directed)   

  result<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                      directed=directed,pValueAdj="BH",
                      pValueCutoff=linkerPValThreshold,
                      communityMethod=communityDetectionMethod,
                      keepIsolatedNodes=keepIsolatedNodes)

  #names(result)

###
  linkerDF<-result$neighborData
  linkerDF<-linkerDF[order(linkerDF$pValueRaw),]
  #linkerDF[linkerDF$pValueFDR<linkerPValThreshold,]
  linkerDF$CGC<-rep("NA",nrow(linkerDF))
  linkerDF$epiEnzyme<-rep("NA",nrow(linkerDF))

  linkerDF[linkerDF$name %in% tmpData$Gene.Symbol,]$CGC<-"YES"
  linkerDF[linkerDF$name %in% tmpData2$HGNC_symbol,]$epiEnzyme<-"YES"
  linkerDF<-linkerDF[linkerDF$pValueFDR<linkerPValThreshold,]

  numOfLinker<-nrow(linkerDF)
  
###

  folderName<-paste("network_numOfLinker",numOfLinker,"communityMethod",communityDetectionMethod,sep="_")
  filePath<-file.path(dataDir,cancerType,"network",folderName)
  if( !file.exists(filePath) ){
    dir.create(filePath,recursive=TRUE)
  }

  write.table(result$netboxOutput[1:3],file=file.path(filePath,"network.sif"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  write.table(result$netboxOutput,file=file.path(filePath,"networkAttr.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
  write.table(result$neighborData,file=file.path(filePath,"neighborList.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

  if( communityDetectionMethod != "none" ){
    fileName<-paste("membership.",communityDetectionMethod,".txt",sep="")
    write.table(result$moduleMembership,file=file.path(filePath,"memb.ebc.txt"),sep="\t",quote=FALSE,col.names=FALSE)
  }
  
  write.table(result$nodeType,file=file.path(filePath,"nodeType.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  
  
  write.table(alterationFreq,file=file.path(filePath,"alterationFreq.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  write.table(alterationList,file=file.path(filePath,"alterationAnno.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

  
  if( numOfLinker>0 ){
    write.table(linkerDF,file=file.path(filePath,"linkerAnno.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
  }
  
  cat(sprintf("NetBox2 output is located at %s\n",filePath))

  netBox2Result<-list(netboxGraph=result$netboxGraph,graphReduced=graphReduced,
                      moduleMembership=result$moduleMembership,
                      sifOutput=result$netboxOutput[1:3],
                      neighborData=result$neighborData)
  
  return(netBox2Result)

}
         