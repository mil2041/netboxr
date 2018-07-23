#' merge copy number call, copy number value, gene expression data frame 
#' 
#' @param CNValueMat GISTIC2 absolute copy number value data frame
#' @param CNCallMat GISTIC2 copy number call data frame
#' @param GEMat gene expression data frame
#' @param geneList gene name of interest
#' @param verbose output detailed information [ TRUE, FALSE ] 
#' 
#' @return a list of data frames 
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' #worDir<-getwd()
#' 
#' #mergedData<-intersectCNGEMat(CNValueMat,CNCallMat,GEMat,geneList)
#' 
#' @concept netboxr
#' @export
intersectCNGEMat<-function(CNValueMat,CNCallMat,GEMat,geneList=NULL,verbose=TRUE){
  
  if(length(nrow(CNValueMat))>0){
    log2CNMat<-log(CNValueMat,base=2)
    log2CNSamples<-colnames(log2CNMat)
    CNValueSamples<-colnames(CNValueMat)
  }else{
    log2CNMat<-{}
    log2CNSample<-{}
    CNValueSamples<-{}
  }
  
  CNCallSamples<-colnames(CNCallMat)
  GESamples<-colnames(GEMat)

  ### take sampleIDs intersection among CNCall and GE 
  commonSamplesID<-intersect(substr(CNValueSamples,1,12),substr(CNCallSamples,1,12))
  commonSamplesID<-intersect(commonSamplesID,substr(GESamples,1,12))
  
  CNValueMat<-CNValueMat[,which(substr(colnames(CNValueMat),1,12) %in% commonSamplesID)]
  CNCallMat<-CNCallMat[,which(substr(colnames(CNCallMat),1,12) %in% commonSamplesID)]
  GEMat<-GEMat[,which(substr(colnames(GEMat),1,12) %in% commonSamplesID)]
  
  if(length(geneList)>0){ 
    log2CNcommonGeneSymbol<-intersect(geneList,rownames(CNValueMat))
    CNCallcommonGeneSymbol<-intersect(geneList,rownames(CNCallMat))
    commonGeneSymbol<-intersect(log2CNcommonGeneSymbol,CNCallcommonGeneSymbol)
    commonGeneSymbol<-intersect(commonGeneSymbol,rownames(GEMat))
    missedGene<-setdiff(geneList,commonGeneSymbol)     
  }else{
    commonGeneSymbol<-intersect(rownames(CNCallMat),rownames(GEMat))    
  }
  
  CNValueMatSub<-CNValueMat[which(rownames(CNValueMat) %in% commonGeneSymbol),]
  CNCallMatSub<-CNCallMat[which(rownames(CNCallMat) %in% commonGeneSymbol),]
  GEMatSub<-GEMat[which(rownames(GEMat) %in% commonGeneSymbol),]
  
  cat(sprintf("Intersection of copy number and RNASeq matrix\n")) 
  cat(sprintf("contains %s samples and %s genes\n",ncol(GEMatSub),nrow(GEMatSub)))
  
  log2CNValueMatSub<-as.data.frame(log(CNValueMatSub,base=2))
  log2CNValueMatSub[log2CNValueMatSub == "-Inf" ]<-NA
  
  log2GEMatSub<-as.data.frame(log(GEMatSub,base=2))
  log2GEMatSub[log2GEMatSub == "-Inf" ]<-NA
  #log2GEMatSub[log2GEMatSub == "-Inf" ]<-(-1)
  
  #log2GEMatSub[log2GEMatSub == "-Inf" ]<-0
  
  mergedData<-list(CNValueMat=CNValueMatSub[order(rownames(CNValueMatSub)),],
                   CNCallMat=CNCallMatSub[order(rownames(CNCallMatSub)),],
                   GEMat=GEMatSub[order(rownames(GEMatSub)),],
                   log2CNValueMat=log2CNValueMatSub[order(rownames(log2CNValueMatSub)),],
                   log2GEMat=log2GEMatSub[order(rownames(log2GEMatSub)),],
                   missedGene=missedGene)
  
  return(mergedData)
  
}
