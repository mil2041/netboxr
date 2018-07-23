#' Building copy number frequency matrix
#' 
#' @param data copy number callk matrix from GISTIC2 
#' @param scale it is (-2,2) in default GISTIC2 convention
#' @param breaks five breaks
#' 
#' @return a data frame of copy number frequency
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
#' #copyNumberFreq<-buildCopyNumberFreqMatrix(wholeCopyNumberMat,
#' #              scale=c(-2,2),breaks=5)/ncol(wholeCopyNumberMat)
#' 
#' @concept netboxr
#' @export
buildCopyNumberFreqMatrix<-function (data,scale,breaks){
  
  nrows<-dim(data)[1]
  ncols<-dim(data)[2]
  scaleMin<-scale[1]
  scaleMax<-scale[2]
  binSize<-( scaleMax - scaleMin + 1 )/breaks
  
  binSeq<-{}
  for(i in 1:breaks) {
    if( i == 1){
      binSeq[i]<-scaleMin
    }else {
      binSeq[i]<-( binSeq[i-1] + binSize )
    }
    
  }
  
  FreqMatrix<-{}
  FreqBinVector<-{}
  for( k in 1: breaks ){
    if( k != breaks ) {         
      FreqBinVector<-( rowSums(data >= binSeq[k]) - rowSums( data >= binSeq[k+1]) )        
    }else {
      FreqBinVector<-rowSums(data >= binSeq[k])
    }
    FreqMatrix<-cbind(FreqMatrix,FreqBinVector)
  }
  
  return(FreqMatrix)
  
}