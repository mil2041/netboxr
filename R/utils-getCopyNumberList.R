#' Retrieve significant copy number list from TCGA study
#' 
#' @param date Default is "latest" 
#' @param cancerType Abbreviation in the TCGA study ["KIRC","GBM"]
#' @param selectedSampleId Subset of TCGA sample id
#' @param workDir Where to save the processed result
#' @param copyNumberFreqThreshold Default is 0.05
#' @param corrAverage Default is "median"
#' @param minRankTestSize Default is 4
#' @param indepThreshold Default is 3
#' @param corrThreshold Default is 0.05
#' @param logFCThreshold Default is 0.5
#' @param verbose Default is TRUE
#' 
#' @return a data frame of data summary
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' selectedSampleId<-NULL
#' #worDir<-getwd()
#' mutSig2CVthreshold<-0.1
#' rareMutationUpperLimit<-0.3
#' rareMutationLowerLimit<-0.1
#' rareMutationFreq<-0.02
#' 
#' #copyNumberList<-getCopyNumberList<-function(date,cancerType, 
#'  #                      selectedSampleId=NA, workDir,
#'  #                      copyNumberFreqThreshold=0.05,
#'  #                      corrAverage="median",
#'  #                      minRankTestSize=4, indepThreshold=3,
#'  #                      corrThreshold=0.05, logFCThreshold=0.5,
#'  #                      verbose=TRUE)
#' 
#' @concept netboxr
#' @export 
#' @importFrom HGNChelper checkGeneSymbols
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
getCopyNumberList<-function(date,cancerType, selectedSampleId=NA, workDir,
                            copyNumberFreqThreshold=0.05,
                            corrAverage="median",
                            minRankTestSize=4, indepThreshold=3,
                            corrThreshold=0.05, logFCThreshold=0.5,
                            verbose=TRUE){

#runDates<-getRunDates()
#date<-runDates[1]
#date<-"2014_10_17"

  cat(sprintf("Processing cancertype: %s\n",cancerType))

  tmpDir<-file.path(dataDir,cancerType,"cna")
  if( !file.exists(tmpDir) ){
    dir.create(tmpDir,recursive=TRUE)
  }


#####

  gistic2CNCallMat<-getGISTIC2Mat(date=date,cancerType=cancerType,workDir=tmpDir)
  dim(gistic2CNCallMat)

  gistic2CNCallMat<-processCNCallMat(gistic2CNCallMat)
  gistic2CNCallMatTumor<-gistic2CNCallMat$tumorMat
  
  if(length(selectedSampleId)==0){
    selectedSampleId<-colnames(gistic2CNCallMatTumor)
  }

  CNCallMat<-gistic2CNCallMatTumor[,which(substr(colnames(gistic2CNCallMatTumor),1,12) %in% substr(selectedSampleId,1,12))]
  
  filePath<-file.path(dataDir,cancerType,"cna")
  save(CNCallMat,file=paste(filePath,"CNCallMat.Rd",sep="/"))

#####

  wholeCopyNumberMat<-CNCallMat
  copyNumberFreq<-buildCopyNumberFreqMatrix(wholeCopyNumberMat,scale=c(-2,2),breaks=5)/ncol(wholeCopyNumberMat)
  copyNumberFreq<-as.data.frame(copyNumberFreq)
  colnames(copyNumberFreq)<-c("Homdel","Hetloss","Diploid","Gain","Amp")

  sampleSize<-ncol(wholeCopyNumberMat)
  #sampleSize

  #filePath<-file.path(PROJHOME,"annotation/UCSC_gene_model")
  #fileName<-"hg19_annotation_by_gene_plot.txt"
  #fileName<-file.path(filePath,fileName)
  
  fileName<-system.file("extdata","ucsc_hg19_annotation_by_gene_plot.txt",package="netboxr")
  
  gene_pos<-read.table(fileName,header=FALSE, sep="\t")
  colnames(gene_pos)<-c("entrezID","geneSymbol","chr","strand","startPos")
  copyNumberFreq$geneSymbol<-rownames(copyNumberFreq)
  dataWithPos<-merge(copyNumberFreq,gene_pos,by="geneSymbol")

  #dataWithPos<-ab
  chrOrder<-c((1:22),"X","Y","M")
  dataWithPos$chr<-factor(dataWithPos$chr,chrOrder,ordered=TRUE)

  dataWithPos<-dataWithPos[order(dataWithPos$chr,dataWithPos$startPos),] 


  txdb<-TxDb.Hsapiens.UCSC.hg19.knownGene
  chromSize<-seqlengths(txdb)[1:25]

  fileName<-paste(cancerType,"_CopyNumber_Freq.png",sep="")
  fileName<-file.path(tmpDir,fileName)
  png(file=fileName,width=1500,height=900)
  plotCopyNumberFreq(dataWithPos=dataWithPos,threshold=copyNumberFreqThreshold,chromSize=chromSize)
  dev.off()

  # set 5% frequency threshold for CNV events
  threshold<-copyNumberFreqThreshold
  keep<-(dataWithPos$Amp > threshold)
  ampSelected<-dataWithPos[keep, ]
  ampGeneNum<-dim(ampSelected)[1]
  table(ampSelected$chr)

  #threshold<-0.05
  keep<-(dataWithPos$Homdel > threshold)
  homoDelSelected<-dataWithPos[keep, ]
  homoDelGeneNum<-dim(homoDelSelected)[1]
  table(homoDelSelected$chr)

  #copyNumberGeneList<-union(ampSelected$GeneSymbol,homoDelSelected$GeneSymbol)
  ampGeneList<-ampSelected$geneSymbol
  delGeneList<-homoDelSelected$geneSymbol

  copyNumberGeneList<-union(ampGeneList,delGeneList)
  geneList<-copyNumberGeneList
  cat(sprintf("At frequency threshold: %s, there are %s Amp and %s Homdel CNV\n",threshold,ampGeneNum,homoDelGeneNum))

  #### GISTIC2 copy number value

  GISTIC2CNValueMat<-getGISTIC2CNValueMat(date=date,cancerType=cancerType,workDir=tmpDir)
  dim(GISTIC2CNValueMat)

  absoluteCNValueMat<-processCNValueMat(GISTIC2CNValueMat)

  #### mRNASeq

  filePath<-file.path(dataDir,cancerType,"geneExpression")
  rawCountsMat<-getRNASeqV2RawCounts(date=date,cancerType=cancerType,workDir=filePath)
  dim(rawCountsMat)

  rawCountsMat<-processRawCountsMat(rawCountsMat,normalization="none")

  gExpMatTumor<-rawCountsMat$tumorMat
  gExpMatNormal<-rawCountsMat$normalMat

  gExpMatTumorSub<-gExpMatTumor[,which(substr(colnames(gExpMatTumor),1,12) %in% substr(selectedSampleId,1,12))]
   
  if(!is.null(gExpMatNormal)){
    gExpMat<-cbind(gExpMatTumorSub,gExpMatNormal)
  
    gExpMat<-processGECountsMat(gExpMat,normalization="TMM")
  
    GEMat<-gExpMat$tumorMat
    GEMatNormal<-gExpMat$normalMat
  }else{
    gExpMat<-gExpMatTumorSub
    
    gExpMat<-processGECountsMat(gExpMat,normalization="TMM")
    
    GEMat<-gExpMat$tumorMat
    GEMatNormal<-NULL
    
  }
  
  
  
  filePath<-file.path(dataDir,cancerType,"geneExpression")
  if( !file.exists(filePath) ){
    dir.create(filePath,recursive=TRUE)
  }

  save(GEMat,file=file.path(filePath,"GEMat.Rd"))
  save(GEMatNormal,file=file.path(filePath,"GEMatNormal.Rd"))

  #### Pre-processing copy number and gene expression matrix

  #rawMat<-processCNCallMat(GISTIC2Mat)

  CNValueMat<-absoluteCNValueMat$tumorMat
  CNCallMat<-gistic2CNCallMat$tumorMat
  GEMat<-gExpMat$tumorMat

  fileName<-"CNCallMat.Rd"
  fileName<-file.path(dataDir,cancerType,"cna",fileName)
  save(CNCallMat,file=fileName)

######
  mergedData<-intersectCNGEMat(CNValueMat,CNCallMat,GEMat,geneList)
  missedGeneNum<-length(mergedData$missedGene)
  cat(sprintf("there are %s genes ignored from merging\n",missedGeneNum))

  filePath<-file.path(dataDir,cancerType,"cna")
  fileName<-paste(cancerType,"_ignored_gene_from_CNCall_GE_merge.txt",sep="")
  fileName<-file.path(filePath,fileName)
  write.table(mergedData$missedGene,file=fileName,sep="\n",row.names = FALSE,col.names = FALSE,quote=FALSE)

  log2CNMat<-mergedData$log2CNValueMat
  CNCallMat<-mergedData$CNCallMat
  log2GEMat<-mergedData$log2GEMat

#### calculate CN - gene expression correlation among extreme and diploid within tumor samples

  filePath<-file.path(dataDir,cancerType)
  if( !file.exists(filePath) ){
    dir.create(filePath,recursive=TRUE)
  }

  fileName<-paste(cancerType,"_corr.Rd",sep="")
  fileName<-file.path(filePath,"CN_GE_corr_MannWhitney",fileName)

  if( !file.exists(fileName) ){
  
    corrTable<-calculateCNGEcorr(log2CNMat,CNCallMat,log2GEMat,ampGeneList,delGeneList,method="MannWhitney",average=corrAverage,rankTestSize=minRankTestSize,workDir=filePath)
    names(corrTable)
  
    corrTestCombined<-rbind(corrTable$ampTest,corrTable$delTest)
    save(corrTestCombined,file=fileName)
  
  }else{
    cat(sprintf("correlation test Rd file exists, skip calculation\n"))
    load(fileName)
  
  }

  # independent filtering
  # remove diploid median less than log(2,base=2)
  #indepThreshold<-3
  corrTestIndep<-corrTestCombined[corrTestCombined$log2DiploidMean>indepThreshold,]

  corrTestIndep$bonferroni<-corrTestIndep$MannWhitneyPValue*nrow(corrTestIndep)
  corrTestIndep[corrTestIndep$bonferroni>=1,]$bonferroni<-1

  corrTestIndep$fdr<-p.adjust(corrTestIndep$MannWhitneyPValue,method="BH")
  corrTestIndep<-corrTestIndep[order(corrTestIndep$copyNumberType,corrTestIndep$fdr),]

  corrTestIndep<-merge(corrTestIndep,dataWithPos,by="geneSymbol")
  #bb<-merge(corrTestIndep,dataWithPos,by="geneSymbol")
  #bb<-bb[order(bb$copyNumberType,bb$fdr),]
  #bb1<-bb[bb$copyNumberType %in% "Amp",]
  #bb2<-bb[bb$copyNumberType %in% "Homdel",]

  #corrThreshold<-0.05
  #logFCthreshold<-0.5
  #filteredGeneList<-corrTestIndep[corrTestIndep$bonferroni>corrThreshold,]$geneSymbol
  #keepedGeneList<-corrTestIndep[corrTestIndep$bonferroni<=corrThreshold,]
  keepedGeneList<-corrTestIndep[corrTestIndep$fdr<=corrThreshold,]
  keep<- keepedGeneList$log2Fold<(-logFCThreshold) | keepedGeneList$log2Fold>(logFCThreshold) 
  keepedGeneList<-keepedGeneList[keep,]
  keepedGeneList<-keepedGeneList[order(keepedGeneList$copyNumberType,keepedGeneList$fdr),]

  filePath<-file.path(dataDir,cancerType,"CN_GE_corr_MannWhitney")

  fileName<-paste(cancerType,"_CNV_corr_test_indep.txt",sep="")
  fileName<-file.path(filePath,fileName)
  write.table(corrTestIndep,file=fileName,quote=FALSE,row.names=FALSE,sep="\t")

  fileName<-paste(cancerType,"_CNV_corr_test_keeped.txt",sep="")
  fileName<-file.path(filePath,fileName)
  write.table(keepedGeneList,file=fileName,quote=FALSE,row.names=FALSE,sep="\t")

  totalCNVnum<-nrow(corrTestCombined)
  removedCNVnum<-(totalCNVnum - nrow(keepedGeneList))

  filteredRatio<-(removedCNVnum/totalCNVnum)
  cat(sprintf("At bonferroni threshold: %s \n",corrThreshold))
  cat(sprintf("At log2FC threshold: %s \n",logFCThreshold))
  cat(sprintf("%s CNV removed / %s CNV total\n",removedCNVnum,totalCNVnum,filteredRatio))
  cat(sprintf("%s CNV passed filteration\n",totalCNVnum-removedCNVnum))
  cat(sprintf("filtration passed ratio: %.2f\n",1-filteredRatio))

  dataWithPosFiltered<-dataWithPos[which(dataWithPos$geneSymbol %in% keepedGeneList$geneSymbol),]

  fileName<-paste(cancerType,"_CopyNumber_Freq_filtered.png",sep="")
  filePath<-file.path(dataDir,cancerType,"CN_GE_corr_MannWhitney")
  fileName<-file.path(filePath,fileName)
  png(file=fileName,width=1024,height=500)
  plotCopyNumberFreq(dataWithPos=dataWithPosFiltered,threshold=copyNumberFreqThreshold,chromSize=chromSize)
  dev.off()

### volcano plot for CN-GE corrleation

  fileName<-paste(cancerType,"_CopyNumber_corr_volcano.png",sep="")
  filePath<-file.path(dataDir,cancerType,"CN_GE_corr_MannWhitney")
  fileName<-file.path(filePath,fileName)

  png(file=fileName,width=1024,height=780)
  ampGroup<-corrTestIndep[corrTestIndep$copyNumberType=="Amp",]
  delGroup<-corrTestIndep[corrTestIndep$copyNumberType=="Homdel",]
  plot(delGroup$log2Fold,-log(delGroup$bonferroni,base=10),pch=20,col="blue", cex=1.5,
     xlab="log2FoldChange [ extreme / diploid ]",ylab="-log10 ( pValue-adjusted )",
     xlim=c(min(delGroup$log2Fold,ampGroup$log2Fold),max(delGroup$log2Fold,ampGroup$log2Fold)),
     ylim=c(0,max(-log(delGroup$bonferroni,base=10),-log(ampGroup$bonferroni,base=10))))

  title(paste("Volcano plot for CNV correlation test [ ",cancerType," ]",sep=""))
  points(ampGroup$log2Fold,-log(ampGroup$bonferroni,base=10),pch=20,col="red",cex=1.5)
  abline(h=(-log(corrThreshold,base=10)),col="green",lwd=3)
  abline(v=c(-logFCThreshold,logFCThreshold),col="grey",lwd=3)

  pValueLegend<-paste(corrThreshold*100,"% pValue",sep="")
  logFCLegend<-paste("logFC = ",logFCThreshold,sep="")
  legend("topright",pch=20,col=c("red","blue","green","grey"),legend=c("Amp","Homdel",pValueLegend,logFCLegend))
  mtext(paste(totalCNVnum-removedCNVnum," passed / ",totalCNVnum," total, rate: ",signif(1-filteredRatio,digit=2),sep=""))
  dev.off()

  head(dataWithPos)

  #copyNumberGeneListFiltered<-corrTestCombined[corrTestCombined$padj<=corrThreshold,]$geneSymbol
  ampGeneListFiltered<-keepedGeneList[keepedGeneList$copyNumberType=="Amp",]$geneSymbol
  delGeneListFiltered<-keepedGeneList[keepedGeneList$copyNumberType=="Homdel",]$geneSymbol

  ampGeneFreq<-copyNumberFreq[as.character(ampGeneListFiltered),]$Amp
  delGeneFreq<-copyNumberFreq[as.character(delGeneListFiltered),]$Homdel
  
  geneList<-union(ampGeneListFiltered,delGeneListFiltered)

  cat(sprintf("making cna data matrix file\n"))
  copyNumberMat<-makeCopyNumberMatrix(CNCallMat,selectedSampleId,geneList)

  fileName<-paste(cancerType,"_copyNumberMatrix.txt",sep="")
  filePath<-file.path(dataDir,cancerType,"cna")
  fileName<-file.path(filePath,fileName)
  write.table(copyNumberMat,file=fileName,sep="\t",row.names = TRUE,col.names = TRUE,quote=FALSE)

  copyNumberList<-list(ampGeneList=ampGeneListFiltered,
                       ampGeneFreq=ampGeneFreq,
                       delGeneList=delGeneListFiltered,
                       delGeneFreq=delGeneFreq,
                       copyNumberMat=copyNumberMat)
  
  return(copyNumberList)

}
