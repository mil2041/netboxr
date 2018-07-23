#' Make copy number frequency plot
#' 
#' @param dataWithPos Default is "latest" 
#' @param threshold Abbreviation in the TCGA study ["KIRC","GBM"]
#' @param chromSize MutSig2CV result data frame
#' @param verbose output detailed information
#' 
#' @return null
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
#' #plotCopyNumberFreq(dataWithPos=dataWithPos,
#' #          threshold=copyNumberFreqThreshold,chromSize=chromSize)
#' 
#' @concept netboxr
#' @export
plotCopyNumberFreq<-function(dataWithPos,threshold,chromSize,verbose=TRUE){
  
  keep<-(dataWithPos$Amp > threshold)
  ampSelected<-dataWithPos[keep, ]
  ampGeneNum<-dim(ampSelected)[1]
  table(ampSelected$chr)
  
  #threshold<-0.05
  keep<-(dataWithPos$Homdel > threshold)
  homoDelSelected<-dataWithPos[keep, ]
  homoDelGeneNum<-dim(homoDelSelected)[1]
  table(homoDelSelected$chr)
  
  # set position, ticks and labels
  dataWithPos$pos<-NA
  
  # re-index, if one chromosome is missing 
  dataWithPos$index<-NA
  ind<-0
  for (i in unique(dataWithPos$chr)){
    ind = ind + 1
    dataWithPos[dataWithPos$chr==i,]$index = ind
  }
  
  ### prepare tickers
  nchr=length(unique(dataWithPos$chr)) 
  
  if (nchr==1) {
    dataWithPos$pos=dataWithPos$startPos
    ticks=floor(length(dataWithPos$startPos))/2+1
    xlabel = paste('chromosome',unique(dataWithPos$chr),'position')
    labs = ticks
  } else {
    ticks = rep(NA,length(unique(dataWithPos$chr))+1)
    ticks[1] = 0
    
    for (i in 1:max(dataWithPos$index)) {
      #dataWithPos[dataWithPos$index==i, ]$pos <- (dataWithPos[dataWithPos$index==i, ]$startPos - dataWithPos[dataWithPos$index==i,]$startPos[1]) +1 +ticks[i]
      #ticks[i+1] <- max(dataWithPos[dataWithPos$index==i,]$pos)
      #cat(sprintf("index:%s\n",i))
      
      dataWithPos[dataWithPos$index==i, ]$pos <- (dataWithPos[dataWithPos$index==i, ]$startPos + ticks[i])
      ticks[i+1]<-sum(as.numeric(chromSize[1:i]))
    }
    
    
    xlabel = 'chromosome'
    labs = c(unique(dataWithPos$chr),'')
  }
  
  # Initialize plot
  
  xmax<-max(dataWithPos$pos) + (1e8)
  xmin<-min(dataWithPos$pos) - (1e8)
  
  ymax<-min(max(dataWithPos[,6]*100) +20, 100)
  ymin<-max(dataWithPos[,2]*100)*(-1) -3

  #ymax<-100
  #ymin<-(-100)
    
  cex.axis=0.8
  
  #threshold<-0.02
  
  #plot(dataWithPos$pos,dataWithPos[,6]*100,type="h",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
  #     xaxt='n',xaxs='i',las=1, yaxt='n',yaxs='i',
  #     ,cex.axis=cex.axis,xlab=xlabel,ylab="Frequency (%)")
  
  #title(main=paste("Copy number alteration: threshold (",threshold*100,"%)",sep=""))
  #mtext(paste("Amplification genes: ",ampGeneNum," homoDeletion genes: ",homoDelGeneNum,sep=""))
  
  #lines(dataWithPos$pos,dataWithPos[,5]*100,type="h",col="grey",ylim=c(-100,100))
  #lines(dataWithPos$pos,-1*(dataWithPos[,3]*100),type="h",col="grey",ylim=c(-100,100))
  #lines(dataWithPos$pos,-1*(dataWithPos[,2]*100),type="h",col="blue",ylim=c(ymin,ymax))
  
  
  #####
  plot(dataWithPos$pos,dataWithPos[,6]*100,type="h",col="red",xlim=c(xmin,xmax),ylim=c(ymin,ymax),
       xaxt='n',xaxs='i',las=1, yaxt='n',yaxs='i',
       ,cex.axis=cex.axis,xlab=xlabel,ylab="Frequency (%)")
  
  title(main=paste("Copy number alteration: threshold (",threshold*100,"%)",sep=""))
  mtext(paste("Amplification genes: ",ampGeneNum," homoDeletion genes: ",homoDelGeneNum,sep=""))
  
  #lines(dataWithPos$pos,dataWithPos[,5]*100,type="h",col="grey",ylim=c(-100,100))
  #lines(dataWithPos$pos,-1*(dataWithPos[,3]*100),type="h",col="grey",ylim=c(-100,100))
  lines(dataWithPos$pos,-1*(dataWithPos[,2]*100),type="h",col="blue",ylim=c(ymin,ymax))
  
  
  #####
  
  abline(v=ticks,col="grey",lty=3)
  abline(h=c(-threshold*100,threshold*100),col="green",lwd=2)
  
  legend("topright", legend=c("Amplification (+2)","Deletion (-2)"), col=c("red","blue"), pch=20,cex=0.9, xpd=TRUE)
  
  # stagger labels
  blank = rep('',length(labs))
  lowerlabs = rep('',length(labs)-1)
  upperlabs = rep('',length(labs)-1)
  
  for (i in 1:length(labs)-1){
    if (i %% 2 == 0){
      lowerlabs[i] = labs[i]
    } else{
      upperlabs[i] = labs[i]
    }
  }
  
  
  if(length(upperlabs)==23){
    upperlabs[23]<-"X"
  }
  
  if(length(lowerlabs)==24){
    lowerlabs[24]<-"Y"
  } 
    
  middleTicks<-0
  for(i in 1:(length(ticks)-1)){
    middleTicks[i]<-(ticks[i] + ticks[i+1] )/2
  }
  
  
  # major x-ticks
  axis(side=1,at=ticks,labels=blank,lwd=0,lwd.ticks=1,cex.axis=cex.axis)
  
  # minor x-ticks
  axis(side=1,at=middleTicks,labels=upperlabs,lwd=0,lwd.ticks=0,cex.axis=cex.axis,line=-0.35)
  axis(side=1,at=middleTicks,labels=lowerlabs,lwd=0,lwd.ticks=0,cex.axis=cex.axis,line=0.35)
  
  
  yTicks<-seq(from=round(ymin,digit=-1),to=round(ymax,digit=-1),by=10)
  axis(side=2,at=yTicks,tck=-0.015,labels=NA)
  axis(side=2,at=yTicks,labels=yTicks,lwd=0, line=0.4,las=1)
  
  minorTicks<-0
  minorTicksNum<-5
  counter<-1
  for(i in 1:(length(yTicks)-1)){
    for(j in 1:(minorTicksNum-1)){
      minorTicks[counter]<-yTicks[i] + ( (yTicks[i+1] -yTicks[i]) / minorTicksNum)*j    
      #cat(counter,j,i,"\n")
      counter<-counter+1
    }
    
  }
  axis(side=2,at=minorTicks,tck=0.010,labels=NA)
  
}