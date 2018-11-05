
filePath<-"~/work/Ekta_lab/netboxr_manuscript/debug"
communityDetectionMethod<-"lec"

linkerDF<-result$neighborData
linkerDF<-linkerDF[order(linkerDF$pValueRaw),]
#linkerDF[linkerDF$pValueFDR<linkerPValThreshold,]
#linkerDF$CGC<-rep("NA",nrow(linkerDF))
#linkerDF$epiEnzyme<-rep("NA",nrow(linkerDF))

#linkerDF[linkerDF$name %in% tmpData$Gene.Symbol,]$CGC<-"YES"
#linkerDF[linkerDF$name %in% tmpData2$HGNC_symbol,]$epiEnzyme<-"YES"
#linkerDF<-linkerDF[linkerDF$pValueFDR<linkerPValThreshold,]

numOfLinker<-nrow(linkerDF)
  
folderName<-paste("network_numOfLinker",numOfLinker,"communityMethod",communityDetectionMethod,sep="_")
#filePath<-file.path(dataDir,cancerType,"network",folderName)
filePath<-file.path(filePath,folderName)
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
