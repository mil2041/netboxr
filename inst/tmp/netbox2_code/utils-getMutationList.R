#' Retrieve significant mutation list from TCGA study
#' 
#' @param date Default is "latest" 
#' @param cancerType Abbreviation in the TCGA study ["KIRC","GBM"]
#' @param selectedSampleId Subset of TCGA sample id
#' @param workDir Where to save the processed result
#' @param mutSig2CVthreshold Default is 0.1
#' @param rareMutationUpperLimit Default is 0.3
#' @param rareMutationLowerLimit Default is 0.1
#' @param rateMutationFreq Default is 0.02
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
getMutationList<-function(date,cancerType, selectedSampleId=NA, workDir,
                          mutSig2CVthreshold=0.1, rareMutationUpperLimit=0.3,
                          rareMutationLowerLimit=0.1, rareMutationFreq=0.02,
                          verbose=TRUE){

  cat(sprintf("Processing cancertype: %s\n",cancerType))

  tmpDir<-file.path(dataDir,cancerType,"mutation")
  if( !file.exists(tmpDir) ){
    dir.create(tmpDir,recursive=TRUE)
  }

#####

  #clinicalDat<-read.table(sampleIdFile,sep="\t",header=TRUE,
  #                      stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')

  #selectedSampleId<-clinicalDat$PATIENT_ID

#####

  filePath<-file.path(dataDir,cancerType,"mutation")
  
  mutationMAF<-getMutationMAF(date,cancerType,workDir = filePath)
  mutationMAF$patientID<-mutationMAF$Tumor_Sample_Barcode
  mutationMAF<-processMutationMAF(mutationMAF)
  
  if(length(selectedSampleId)==0){
    selectedSampleId<-unique(substr(mutationMAF$patientID,1,12))
  }
  
#####  
  
  mutSig2CVMat<-getMutSig2CVMat(date=date,cancerType=cancerType,workDir=tmpDir)
  dim(mutSig2CVMat)

  sampleSize<-length(selectedSampleId)
  cat(sprintf("There are %s samples with mutation data\n",sampleSize))

  filePath<-file.path(dataDir,cancerType,"mutation")
  fileName<-paste(cancerType,"_Mutation_Freq_01.png",sep="")  
  fileName<-file.path(filePath,fileName)
  png(file=fileName,width=1024,height=500)
  plotMutationFreq(cancerType = cancerType,mat=mutSig2CVMat,sampleSize=sampleSize,threshold=0.1)
  dev.off()

  filePath<-file.path(dataDir,cancerType,"mutation")
  fileName<-paste(cancerType,"_Mutation_Freq_02.png",sep="")  
  fileName<-file.path(filePath,fileName)
  png(file=fileName,width=1024,height=500)
  plotMutationFreq(cancerType = cancerType,mat=mutSig2CVMat,sampleSize=sampleSize,threshold=0.2)
  dev.off()

  filePath<-file.path(dataDir,cancerType,"mutation")
  fileName<-paste(cancerType,"_Mutation_Freq_03.png",sep="")  
  fileName<-file.path(filePath,fileName)
  png(file=fileName,width=1024,height=500)
  plotMutationFreq(cancerType = cancerType,mat=mutSig2CVMat,sampleSize=sampleSize,threshold=0.3)
  dev.off()

  mutSig2CVMat$sampleSize<-rep(sampleSize,nrow(mutSig2CVMat))
  mutSig2CVMat$freq<-(mutSig2CVMat$npat/mutSig2CVMat$sampleSize)

  checkTable<-suppressWarnings(checkGeneSymbols(as.character(mutSig2CVMat$gene)))
  colnames(checkTable)<-c("gene","Approved","suggestedSymbol")

  mutSig2CVMat<-merge(mutSig2CVMat,checkTable,by="gene")
  mutSig2CVMat<-mutSig2CVMat[!is.na(mutSig2CVMat$suggestedSymbol),]

## check 573 cancer census genes
  #filePath<-file.path(PROJHOME,"annotation/Cosmic_census_list")
  #fileName<-"cosmic_CGC_0128_2016.txt"
  #fileName<-file.path(filePath,fileName)
  
  fileName<-system.file("extdata","cosmic_CGC_0128_2016.txt",package="netboxr")
  
  tmpData<-read.table(fileName,sep="\t",header=TRUE,
                    stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')

  dd<-mutSig2CVMat[ mutSig2CVMat$suggestedSymbol %in% tmpData$Gene.Symbol,]
  #nn<-dd[order(-dd$freq),]
  #head(nn[nn$freq>=0.02,])
  dd1<-dd[dd$q<=rareMutationUpperLimit & dd$q>rareMutationLowerLimit,]
  dd2<-dd[dd$q>rareMutationUpperLimit & dd$freq>=rareMutationFreq, ]
  MutationExtraCGC<-rbind(dd1,dd2)

## check 720 epi enzymes
  #filePath<-file.path(PROJHOME,"annotation/epigenetic_gene_list/epifactors")
  #fileName<-"EpiGenes_main_genes.txt"
  #fileName<-file.path(filePath,fileName)
  
  fileName<-system.file("extdata","EpiGenes_main_genes.txt",package="netboxr")
  
  tmpData2<-read.table(fileName,sep="\t",header=TRUE,
                     stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')

  dd<-mutSig2CVMat[ mutSig2CVMat$gene %in% tmpData2$HGNC_symbol,]
  dd1<-dd[dd$q<=rareMutationUpperLimit & dd$q>rareMutationLowerLimit,]
  dd2<-dd[dd$q>rareMutationUpperLimit & dd$freq>=rareMutationFreq, ]
  MutationExtraEpi<-rbind(dd1,dd2)

####
  if(FALSE){
  
    #filePath<-file.path(PROJHOME,"annotation/kinase_list")
    #fileName<-"kinbase_03102015.txt"
    #fileName<-file.path(filePath,fileName)
  
    fileName<-system.file("extdata","kinbase_03102015.txt",package="netboxr")
    
    tmpData3<-read.table(fileName,sep="\t",header=TRUE,
                       stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')
  
    tmpData3<-tmpData3[,-1]
  
    checkTable<-checkGeneSymbols(as.character(tmpData3$Gene))
    colnames(checkTable)<-c("Gene","Approved","suggestedSymbol")
  
    tmpData3<-merge(tmpData3,checkTable,by="Gene")
    tmpData3<-tmpData3[!is.na(tmpData3$suggestedSymbol),]
  
    dd<-mutSig2CVMat[ mutSig2CVMat$suggestedSymbol %in% tmpData3$suggestedSymbol,]
    dd1<-dd[dd$q<=rareMutationUpperLimit & dd$q>rareMutationLowerLimit,]
    dd2<-dd[dd$q>rareMutationUpperLimit & dd$freq>=rareMutationFreq, ]
    MutationExtraKinase<-rbind(dd1,dd2)
  
  }


  MutationExtra<-unique(rbind(MutationExtraCGC,MutationExtraEpi))
  MutationExtra<-MutationExtra[order(MutationExtra$q,-MutationExtra$freq),]
  filePath<-file.path(dataDir,cancerType,"mutation")
  fileName<-paste(cancerType,"_extra_mutation.txt",sep="")
  fileName<-file.path(filePath,fileName)
  write.table(MutationExtra,file=fileName,sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)

  #intersect(mutSig2CVMat$gene, tmpData$geneSymbol)
  threshold<-mutSig2CVthreshold
  MutationThreshold<-threshold
  keep<-mutSig2CVMat$q<MutationThreshold
  MutationSig<-mutSig2CVMat[keep,]
  dim(MutationSig)

  filePath<-file.path(dataDir,cancerType,"mutation")
  fileName<-paste(cancerType,"_MutSigCV2_default_mutation.txt",sep="")
  fileName<-file.path(filePath,fileName)
  write.table(MutationSig,file=fileName,sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)


  selectedMutation<-rbind(MutationSig,MutationExtra)
  selectedMutation<-selectedMutation[order(selectedMutation$q,-selectedMutation$freq),]
  save(selectedMutation,file=file.path(filePath,"selectedMutation.Rd"))
  #load(paste(workDir,"selectedMutation.Rd",sep="/"))

  if("MLL2" %in% selectedMutation$gene){
    selectedMutation[selectedMutation$gene %in% "MLL2",]$suggestedSymbol<-"KMT2B"
  }
  
  #sigGenesMutation<-MutationSig
  sigGenesMutation<-selectedMutation

  filePath<-file.path(dataDir,cancerType,"mutation")
  fileName<-paste(cancerType,"_all_selected_mutation.txt",sep="")
  fileName<-file.path(filePath,fileName)
  write.table(selectedMutation,file=fileName,sep="\t",row.names = FALSE,col.names = TRUE,quote=FALSE)

#####

  cat(sprintf("making mutation data matrix file\n"))
  geneList<-sigGenesMutation$gene
  
  mutationMat<-makeMutationMatrix(mutationMAF,selectedSampleId,geneList)
  
  filePath<-file.path(dataDir,cancerType,"mutation")
  fileName<-paste(cancerType,"_mutationMatrix.txt",sep="")
  fileName<-file.path(filePath,fileName)
  write.table(mutationMat,file=fileName,sep="\t",row.names = TRUE,col.names = TRUE,quote=FALSE)

  
  mutationResult<-list(mutationList=sigGenesMutation$gene,
                       mutationFreq=sigGenesMutation$freq,
                       mutationMat=mutationMat)
  
  return(mutationResult)

}
