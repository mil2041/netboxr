#' Make MutSig2CV result plot
#' 
#' @param date Default is "latest" 
#' @param cancerType Abbreviation in the TCGA study ["KIRC","GBM"]
#' @param mat MutSig2CV result data frame
#' @param sampleSize 
#' @param threshold Default is 0.05
#' @param workDir Where to save the processed result
#' @param verbose output detailed information
#' 
#' @return a data frame of annotated mutation MAF file
#'   
#' @examples
#' cancerType<-"KIRC"
#' date<-"2015_08_21" 
#' #mutSig2CVMat<-getMutSig2CVMat(date=date,cancerType=cancerType,workDir=tmpDir)
#' #dim(mutSig2CVMat)
#' sampleSize<-333
#' threshold<-0.05
#' #worDir<-getwd()
#' 
#' #mutationMAF<-plotMutationFreq(cancerType,mat,sampleSize,threshold,workDir)
#' 
#' @concept netboxr
#' @export
plotMutationFreq<-function(cancerType,mat,sampleSize,threshold,workDir,verbose){

  MutSig2CVMat<-mat

  MutSig2CVMat$sampleSize<-rep(sampleSize,nrow(MutSig2CVMat))
  MutSig2CVMat$freq<-(MutSig2CVMat$npat/MutSig2CVMat$sampleSize)

  MutationThreshold<-threshold
  keep<-MutSig2CVMat$q<MutationThreshold
  MutationSig<-MutSig2CVMat[keep,]
  dim(MutationSig)

  cat(sprintf("At q-Value: %s, there are %s genes from %s patients\n",
            MutationThreshold,nrow(MutationSig),MutSig2CVMat$sampleSize[1]))

  MutationVector<-MutationSig$freq
  names(MutationVector)<-MutationSig$gene

  x<-barplot(MutationVector*100,col="grey",xaxt="n",ylab="Frequency(%)")
  title(paste(cancerType," Mutation Plot, MutSig2CV q-Value:",
              MutationThreshold,sep=""))
  mtext(paste(nrow(MutationSig)," genes from ",sampleSize," patients",sep=""))
  labs<-as.character(MutationSig$gene)
  #text(x=x, y=-1.25,labels=labs, pos=1,xpd=TRUE,srt=45)
  #axis(side=1,at=x,labels=labs,cex.axis=1,srt=45)
  axis(side=1,at=x,labels=labs,lwd=0,lwd.ticks=0,cex.axis=0.8,line=-0.35,las=2)
  abline(h=5,col="green",lwd=3)

  legend("topright", 
       legend = c("5 % frequency"), 
       fill = c("green"))

}



