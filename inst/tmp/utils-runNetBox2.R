library(HGNChelper)
library(gplots)
library(netboxr)

#####


#####


runNetBox2<-function(mutationList,
                     copyNumberList,
                     methylationList,
                     pathwayCommonDb,
                     directed,
                     linkerPValThreshold,
                     communityDetectionMethod,
                     keepIsolatedNodes,
                     verbose=TRUE){

mutationList<-sigGenesMutation$gene


venn(list(ampCN=ampGeneListFiltered,delCN=delGeneListFiltered))

fileName<-paste(cancerType,"_venn_geneList.png",sep="")
filePath<-file.path(dataDir,cancerType)
fileName<-file.path(filePath,fileName)
png(file=,width=1024,height=780)
#venn(list(mutation=sigGenesMutation$gene,ampCN=ampGeneListFiltered,delCN=delGeneListFiltered))
venn(list(mutation=sigGenesMutation$gene,ampCN=ampGeneListFiltered,delCN=delGeneListFiltered,hyperMet=methylatedList))
dev.off()

intersect(methylatedList,ampGeneListFiltered)
intersect(sigGenesMutation$gene,delGeneListFiltered)


#venn(list(mutation=sigGenesMutation$gene,ampCN=ampGeneListFiltered,delCN=delGeneListFiltered,met=metGeneList))

copyNumberGeneListFiltered<-union(ampGeneListFiltered,delGeneListFiltered)
#copyNumberGeneListFiltered<-setdiff(copyNumberGeneList,log2CNGEnotCorrelatedPearson$GeneSymbol)

#copyNumberList<-copyNumberGeneListFiltered
#geneList<-union(mutationList,copyNumberList)
#length(geneList)


#library(org.Hs.eg.db)
#library(annotate)


checkTable<-suppressWarnings(checkGeneSymbols(as.character(mutationList)))
checkTable$freq<-sigGenesMutation$freq

# convert to latest HGNC geneSymbol and remove obselete geneSymbol
checkTable<-checkTable[!is.na(checkTable[,3]),]
mutationList<-checkTable[,3]
mutationFreq<-checkTable[,4]
mutationDF<-data.frame(mutationList,mutationFreq,rep("Mutation",length(mutationList)))
colnames(mutationDF)<-c("geneSymbol","freq","type")

if(length(ampGeneListFiltered)>0){
checkTable<-suppressWarnings(checkGeneSymbols(as.character(ampGeneListFiltered)))
checkTable$freq<-copyNumberFreq[as.character(ampGeneListFiltered),]$Amp
# convert to latest HGNC geneSymbol and remove obselete geneSymbol
checkTable<-checkTable[!is.na(checkTable[,3]),]
ampGeneListFiltered<-checkTable[,3]
ampGeneFreq<-checkTable[,4]
ampDF<-data.frame(ampGeneListFiltered,ampGeneFreq,rep("Amp",length(ampGeneListFiltered)))
colnames(ampDF)<-c("geneSymbol","freq","type")
}else{ ampDF<-{}}

if(length(delGeneListFiltered)>0){

checkTable<-suppressWarnings(checkGeneSymbols(as.character(delGeneListFiltered)))
checkTable$freq<-copyNumberFreq[as.character(delGeneListFiltered),]$Homdel

# convert to latest HGNC geneSymbol and remove obselete geneSymbol
checkTable<-checkTable[!is.na(checkTable[,3]),]
delGeneListFiltered<-checkTable[,3]
delGeneFreq<-checkTable[,4]
delDF<-data.frame(delGeneListFiltered,delGeneFreq,rep("homDel",length(delGeneListFiltered)))
colnames(delDF)<-c("geneSymbol","freq","type")
}else{ delDF<-{}}


alterationFreq<-rbind(mutationDF,ampDF,delDF)
alterationFreq<-alterationFreq[order(alterationFreq$geneSymbol,-alterationFreq$freq),]
#keep<-!duplicated(alterationFreq$geneSymbol)
#alterationFreq<-alterationFreq[keep,]


hyperMetDF<-data.frame(methylatedList,rep(40,length(methylatedList)),rep("hyperMet",length(methylatedList)))
colnames(hyperMetDF)<-c("geneSymbol","freq","type")

#####

## check 573 cancer census genes
filePath<-file.path(PROJHOME,"annotation/Cosmic_census_list")
fileName<-"cosmic_CGC_0128_2016.txt"
fileName<-file.path(filePath,fileName)

tmpData<-read.table(fileName,sep="\t",header=TRUE,
                    stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')


## check 720 epi enzymes
filePath<-file.path(PROJHOME,"annotation/epigenetic_gene_list/epifactors")
fileName<-"EpiGenes_main_genes.txt"
fileName<-file.path(filePath,fileName)

tmpData2<-read.table(fileName,sep="\t",header=TRUE,
                     stringsAsFactors=FALSE,fill=TRUE,quote=NULL, comment='')


alterationList<-rbind(mutationDF,ampDF,delDF,hyperMetDF)
alterationList$CGC<-rep("NA",nrow(alterationList))
alterationList$epiEnzyme<-rep("NA",nrow(alterationList))

alterationList[alterationList$geneSymbol %in% tmpData$Gene.Symbol,]$CGC<-"YES"
alterationList[alterationList$geneSymbol %in% tmpData2$HGNC_symbol,]$epiEnzyme<-"YES"
  
  
#####  

#fusionList<-c("ERG","ETV1","ETV4","FLI1")

#####
  
copyNumberGeneListFiltered<-union(ampGeneListFiltered,delGeneListFiltered)
#copyNumberGeneListFiltered<-setdiff(copyNumberGeneList,log2CNGEnotCorrelatedPearson$GeneSymbol)

copyNumberList<-copyNumberGeneListFiltered
geneList<-union(mutationList,copyNumberList)
geneList<-union(geneList,methylatedList)
#geneList<-union(geneList,fusionList)
length(geneList)



#setwd("~/work/Sander_lab/TCGA_data/BRCA/TCGA/CN_GE_correlation")
#load("geneList.Rd")


#library(HGNChelper)
data(PC2V7_02212016_sif)
network<-PC2V7_02212016_sif
dim(network)
names(network)
#colnames(network)<-c("geneA","interaction","geneB","databaseName")

data(bioGridNetworkHomoSapiens34133)
networkBioGrid<-bioGridNetworkHomoSapiens34133
dim(networkBioGrid)
names(networkBioGrid)

networkBioGridE2<-networkBioGrid[networkBioGrid$EVIDENCE_NUM>=2,]
dim(networkBioGridE2)

#networkBioGridE3<-networkBioGrid[networkBioGrid$EVIDENCE_NUM==1,]
#trend<-table(networkBioGridE3$YEAR)
#plot(trend)
#networkBioGridE5<-networkBioGridE3[networkBioGridE3$YEAR>=2010,]
#dim(networkBioGridE5)


#networkBioGridE6<-rbind(networkBioGridE2,networkBioGridE5)

selectedNetwork<-rbind(network,networkBioGridE2[,1:5])

fileName<-"netbox2Network_w_Biogrid_e2.sif"
filePath<-file.path(dataDir)
fileName<-file.path(filePath,fileName)
write.table(selectedNetwork,file=fileName,sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)


#selectedNetwork2<-rbind(network,networkBioGridE6[,1:5])

##
# simpilify function confilict in IRange and igraph

# communityMethod=c(ebc,lec)
# pValueAdj=c(BH,Bonferroni)
#result<-geneConnector(geneList=geneList,sifNetwork=network,pValueAdj="Bonferroni",pValueCutoff=0.05,communityMethod="lec",keepIsolatedNodes=TRUE)

#TFList<-c("HIF1A","STAT3","CEBPA","CEBPB","RUNX1","FOSL2","BHLHB2")
sifNetwork<-selectedNetwork
#sifNetwork2<-selectedNetwork2
#sifNetwork<-network
graphReduced<-networkSimplify(sifNetwork,directed = FALSE)   
#graphReduced2<-networkSimplify(sifNetwork2,directed = FALSE)   
#threshold<-0.05
result<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                      directed=FALSE,pValueAdj="BH",
                      pValueCutoff=linkerPValThreshold,
                      communityMethod=communityDetectionMethod,
                      keepIsolatedNodes=keepIsolatedNodes)


#result<-geneConnector(geneList=geneList,networkGraph=graphReduced,directed=FALSE,pValueAdj="BH",
#                       pValueCutoff=threshold,communityMethod="ebc",keepIsolatedNodes=FALSE)




names(result)

###
linkerDF<-result$neighborData
linkerDF<-linkerDF[order(linkerDF$pValueRaw),]

linkerDF[linkerDF$pValueFDR<threshold,]
linkerDF$CGC<-rep("NA",nrow(linkerDF))
linkerDF$epiEnzyme<-rep("NA",nrow(linkerDF))

linkerDF[linkerDF$name %in% tmpData$Gene.Symbol,]$CGC<-"YES"
linkerDF[linkerDF$name %in% tmpData2$HGNC_symbol,]$epiEnzyme<-"YES"
linkerDF<-linkerDF[linkerDF$pValueFDR<threshold,]




###

#localNullModel(result$netboxGraph,iterations=100)

#globeNullModel(netboxGraph =result$netboxGraph ,networkGraph =graphReduced, iterations = 50,numOfGenes = 370)


names(result)

#setwd("~/work/Sander_lab/TCGA_data/BRCA/TCGA/netboxr_Result")
#dataDir<-paste("~/work/Sander_lab/netboxr_result",cancerType,"network/network_wLinker_woHPRD",sep="/")
#if( !file.exists(paste(dataDir,sep="/")) ){
#  dir.create(paste(dataDir,sep=""),recursive=TRUE)
#}

filePath<-file.path(dataDir,cancerType,"network","network_wLinker_wBioGrid_e2_no2_ebc")
if( !file.exists(filePath) ){
  dir.create(filePath,recursive=TRUE)
}


write.table(result$netboxOutput[1:3],file=file.path(filePath,"network.sif"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(result$netboxOutput,file=file.path(filePath,"networkAttr.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(result$neighborData,file=file.path(filePath,"neighborList.txt"),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
write.table(result$moduleMembership,file=file.path(filePath,"memb.ebc.txt"),sep="\t",quote=FALSE,col.names=FALSE)
write.table(result$nodeType,file=file.path(filePath,"nodeType.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(alterationFreq,file=file.path(filePath,"alterationFreq.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(alterationList,file=file.path(filePath,"alterationAnno.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(linkerDF,file=file.path(filePath,"linkerAnno.txt"),sep="\t",quote=FALSE,col.names=FALSE,row.names=FALSE)

cat(sprintf("NetBox2 output is located at %s\n",filePath))

}
         