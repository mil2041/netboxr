data(netbox2010)

geneList<-netbox2010$geneList
sifNetwork<-netbox2010$network

threshold<-0.05

graphReduced<-networkSimplify(sifNetwork,directed = FALSE)
result<-geneConnector(geneList=geneList,networkGraph=graphReduced,directed=FALSE,pValueAdj="BH",
            pValueCutoff=0.05,communityMethod="lec",keepIsolatedNodes=FALSE)

names(result)

linkerDF<-result$neighborData
linkerDF[linkerDF$pValueFDR<threshold,]

localNullModel(result$netboxGraph,iterations=200)

globalNullModel(netboxGraph =result$netboxGraph ,networkGraph =graphReduced, iterations = 30,numOfGenes = 274)




write.table(result$netboxOutput,file="network.sif",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(result$neighborData,file="neighborList.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(result$moduleMembership,file="memb.ebc.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(result$nodeType,file="nodeType.txt",sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

######

extendedSifFileNames<-getPC2NetworkName(fileType="EXTENDED_BINARY_SIF")
selectedFileName<-extendedSifFileNames[grepl("All",extendedSifFileNames)]
extendedSifNetwork <- downloadPc2(selectedFileName=selectedFileName)

names(extendedSifNetwork$edges)

interactionTypes<-unique(extendedSifNetwork$edges$INTERACTION_TYPE)
interactionDataSource<-unique(extendedSifNetwork$edges$INTERACTION_DATA_SOURCE)
