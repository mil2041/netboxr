#' Pre-processing mutation maf file
#' 
#' @param mutationMAF mutation maf file
#' @param verbose Defaukt is TRUE
#' 
#' @return a list of data summary
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' selectedSampleID<-NULL
#' #worDir<-getwd()
#' mutSig2CVthreshold<-0.1
#' rareMutationUpperLimit<-0.3
#' rareMutationLowerLimit<-0.1
#' rareMutationFreq<-0.02
#' 
#' #mutationList<-getMutationList(date,cancerType, selectedSampleId, workDir,
#  #              mutSig2CVthreshold, rareMutationUpperLimit,
#  #              rareMutationLowerLimit, rareMutationFreq)
#' 
#' @concept netboxr
#' @export 
#' @importFrom HGNChelper checkGeneSymbols
processMutationMAF<-function(mutationMAF,verbose=TRUE){

  mutationMAF$patientID<-mutationMAF$Tumor_Sample_Barcode

  sampleIDs<-as.character(unique(mutationMAF$patientID))

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

  # D for DNA and W for whole-genome amplification
  sampleIDs2 <- as.character(samplesDat[,5])
  sampleIDs2 <- substr(sampleIDs2,3,3)

  normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]

  ## the TCGA guide said tumor samples are ID<10, 
  ## but we only use 01A for TP and WXS data
  tumorSamples <- rownames(samplesDat)[sampleIDs1 < 2 & sampleIDs2 =="D" ]

  tumorSamplesCheck<-table(substr(tumorSamples,1,16))

  select<-tumorSamplesCheck==1
  select1<-substr(tumorSamples,1,16) %in% names(tumorSamplesCheck[select])
  tumorSamplesPart1<-tumorSamples[select1]

  select<-tumorSamplesCheck>1
  samplesCheck<-names(tumorSamplesCheck[select])

  tumorSamplesPart2<-0
  for(j in 1:length(samplesCheck)){
  
    select2<-substr(tumorSamples,1,16) %in% samplesCheck[j]
    select2IDs<-tumorSamples[select2]
  
    select2MAFLength<-0
    for (i in 1:length(select2IDs)){
      pick<-mutationMAF$patient %in% select2IDs[i]
      select2MAFLength[i]<-nrow(mutationMAF[pick,])
    }
  
    tumorSamplesPart2[j]<-select2IDs[which.max(select2MAFLength)]
  
  }

  tumorSamples<-sort(append(tumorSamplesPart1,tumorSamplesPart2,after=length(tumorSamplesPart1)))

  sampleSize<-length(tumorSamples)


  if(length(tumorSamples)==0){
    tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
  }

  mutationMAFfixed<-mutationMAF[which(mutationMAF$patientID %in% tumorSamples),]

  return(mutationMAFfixed)

}