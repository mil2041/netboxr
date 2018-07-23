#' Process RNASeqV2 raw count read using edgeR TMM method
#' 
#' @param rawCountsMat raw counts data frame from TCGA
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
#' #rawCountsMat<-processRawCountsMat(rawCountsMat,normalization="TMM")
#' 
#' @concept netboxr
#' @export
#' @import edgeR
processRawCountsMat<-function(rawCountsMat,normalization="TMM"){
  
  cat(sprintf("Pre-processing RNASeqV2 raw counts matrix\n"))
  
  rnaseq_probe_Str<-rawCountsMat[,1][2:nrow(rawCountsMat)]
  out <- strsplit(as.character(rnaseq_probe_Str),"|",fixed=TRUE) 
  rnaseq_probe<-as.data.frame(do.call(rbind,out))
  colnames(rnaseq_probe)<-c("geneSymbol","geneID")
  GEgeneSymbolTable<-rnaseq_probe
  
  GEgeneSymbolTable<-data.frame(lapply(GEgeneSymbolTable, as.character), stringsAsFactors=FALSE)
  GEgeneSymbolTable[GEgeneSymbolTable$geneID %in% "728661",]$geneSymbol<-"SLC35E2B"
  
  totalSampleIDs<-colnames(rawCountsMat)[2:ncol(rawCountsMat)]
  sampleSize<-length(totalSampleIDs)/3
  keep<-seq(1,by=1,len=sampleSize)*3-1
  rawCountsMat<-rawCountsMat[-1,keep]
  
  # bug TCGA-AX-A1C7-01A-11R-A137-07.1 [UCEC]
  keep<-!duplicated(substr(colnames(rawCountsMat),1,28))
  rawCountsMat<-rawCountsMat[,keep]
  
  
  sampleIDs<-unique(colnames(rawCountsMat))
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
  
  
  #keep<-seq(1,by=1,len=length(sampleIDs))*3-1
  
  #remove character lable row
  #rawCountsMat<-rawCountsMat[-1,keep]
  # convert back to numeric column type
  rawCountsMat<-as.data.frame(sapply(rawCountsMat,as.numeric))
  
  criteria1<-duplicated(GEgeneSymbolTable$geneSymbol)
  criteria2<-GEgeneSymbolTable$geneSymbol %in% "?"    
  keep<- !(criteria1 | criteria2)
  keepGeneSymbol<-GEgeneSymbolTable$geneSymbol[keep]
  GEgeneSymbolTable<-GEgeneSymbolTable[keep,]
    
  rawCountsMat<-round(rawCountsMat)
  rawCountsMat<-rawCountsMat[keep,]
  rownames(rawCountsMat)<-keepGeneSymbol
  
  
  
  #rownames(rawCountsMat)<-rnaseq_probe$geneID
  
  
  
  
  #rownames(rawCountsMat)<-rnaseq_probe$entrezID
  
  if( normalization == "TMM"){
  #sampleLibSize<-colSums(rawCountsMat)
  cat(sprintf("Normalizing RNASeq raw counts matrix by edgeR TMM method\n"))
  
  normFactors<-calcNormFactors(rawCountsMat, method="TMM")
  scale<-(colSums(rawCountsMat))*(normFactors)
  normalizedCountsMat<-as.data.frame(round(t(t(rawCountsMat)/scale)*mean(scale)))

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
     GEMatTumor<-rawCountsMat[,GESamplesTumor]
     if( length(normalSamples)>0 ){
         #GESamplesNormal<-gsub("-",".",normalSamples)
         #GEMatNormal<-rawCountsMat[,gsub("\\.","-",GESamplesNormal)]
         GEMatNormal<-rawCountsMat[,GESamplesNormal]
         
     }else{
         GEMatNormal<-NULL
     }
     
  }
  
  cat(sprintf("The file contains %s tumor samples and %s normal samples\n",length(tumorSamples),length(normalSamples)))
  
  GEMat<-list(tumorMat=GEMatTumor,normalMat=GEMatNormal,geneSymbolTable=GEgeneSymbolTable)
  
  return(GEMat)
  
}
