#' Process RNASeqV2 raw count read using edgeR TMM method
#' 
#' @param gExpMat gene expression data frame 
#' @param normalization Abbreviation in the TCGA study ["TMM","none"]
#' @param verbose output detailed information [ TRUE, FALSE ] 
#' 
#' @return a data frame of RNASeq read counts
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' #worDir<-getwd()
#' 
#' #gExpMat<-processGECountsMat(gExpMat,normalization="TMM")
#' 
#' @concept netboxr
#' @export
#' @import edgeR
processGECountsMat<-function(gExpMat,normalization="TMM",verbose=TRUE){
  
  sampleIDs<-unique(colnames(gExpMat))
  length(sampleIDs)
  
  #file_fn<-"RSEM_genes_sample_id.txt"
  #id_df<-read.table(file_fn,sep="\t",header=FALSE)
  
  ### code from RTCGAToolbox.R
  #sampleIDs<-as.character(id_df[,1])
  samplesDat <- data.frame(matrix(nrow=length(sampleIDs),ncol=7))
  rownames(samplesDat) <- sampleIDs
  for(j in 1:length(sampleIDs))
  {
    tmpRow <- unlist(strsplit(sampleIDs[j],split="-"))
    samplesDat[sampleIDs[j],] <- tmpRow
  }
  sampleIDs1 <- as.character(samplesDat[,4])
  sampleIDs1 <- substr(sampleIDs1,1,nchar(sampleIDs1)-1)
  sampleIDs1 <- as.numeric(sampleIDs1)
  normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
  
  ## the TCGA guide said tumor samples are ID<10, but we only use 01A  
  tumorSamples <- rownames(samplesDat)[sampleIDs1 < 2]
  
  if(length(tumorSamples)==0){
    tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
  }
  
  
  #### 
  #GESamplesTumor<-gsub("-",".",tumorSamples)
  
  #if( length(normalSamples)>0 ){
  #  GESamplesNormal<-gsub("-",".",normalSamples)
  #}
  
  GESamplesTumor<-tumorSamples
  
  if( length(normalSamples)>0 ){
    GESamplesNormal<-normalSamples
  }
  
  
  if( normalization == "TMM"){
  #sampleLibSize<-colSums(rawCountsMat)
  cat(sprintf("Normalizing RNASeq raw counts matrix by edgeR TMM method\n"))
  
  normFactors<-calcNormFactors(gExpMat, method="TMM")
  scale<-(colSums(gExpMat))*(normFactors)
  normalizedCountsMat<-as.data.frame(round(t(t(gExpMat)/scale)*mean(scale)))

  #GEMatTumor<-normalizedCountsMat[,gsub("\\.","-",GESamplesTumor)]
  GEMatTumor<-normalizedCountsMat[,GESamplesTumor]
  
    if( length(normalSamples)>0 ){
      #GESamplesNormal<-gsub("-",".",normalSamples)
      #GEMatNormal<-normalizedCountsMat[,gsub("\\.","-",GESamplesNormal)]
      GEMatNormal<-normalizedCountsMat[,GESamplesNormal]
      
      
    }else{
      GEMatNormal<-NULL
    }
   
  }
  
  if( normalization == "none" ){
     cat(sprintf("RNASeq raw counts matrix without normalization\n"))
     #GEMatTumor<-rawCountsMat[,gsub("\\.","-",GESamplesTumor)]
     GEMatTumor<-gExpMat[,GESamplesTumor]
     if( length(normalSamples)>0 ){
         #GESamplesNormal<-gsub("-",".",normalSamples)
         #GEMatNormal<-rawCountsMat[,gsub("\\.","-",GESamplesNormal)]
         GEMatNormal<-gExpMat[,GESamplesNormal]
         
     }else{
         GEMatNormal<-NULL
     }
     
  }
  
  cat(sprintf("The file contains %s tumor samples and %s normal samples\n",length(tumorSamples),length(normalSamples)))
  
  GEMat<-list(tumorMat=GEMatTumor,normalMat=GEMatNormal)
  
  return(GEMat)
  
}
