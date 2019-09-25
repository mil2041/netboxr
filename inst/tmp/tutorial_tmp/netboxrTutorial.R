## ----knitrSetup, include=FALSE---------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", fig.width=12, fig.height=12, tidy=TRUE)

## ----style, include=FALSE, echo=FALSE, results='asis'----------------------
BiocStyle::markdown()

## ----installNetBoxr, eval=FALSE--------------------------------------------
#  library(remotes)
#  install_github(repo="mil2041/netboxr", ref="master", build_vignette=TRUE)

## ----loadLibrary, message=FALSE, warning=FALSE-----------------------------
library(netboxr)

## ----searchHelp, eval=FALSE, tidy=FALSE------------------------------------
#  help.search("netboxr")

## ----showHelp, eval=FALSE, tidy=FALSE--------------------------------------
#  help(geneConnector)
#  ?geneConnector

## ----netboxrExampleNetwork-------------------------------------------------
data(netbox2010)
sifNetwork <- netbox2010$network
graphReduced <- networkSimplify(sifNetwork,directed = FALSE)      

## ----netboxrExampleGene----------------------------------------------------
geneList <- as.character(netbox2010$geneList) 
length(geneList)

## ----netboxrExampleGeneConnector, fig.width=12, fig.height=12--------------
## Use Benjamini-Hochberg method to do multiple hypothesis 
## correction for linker candidates.

## Use edge-betweeness method to detect community structure in the network. 
threshold <- 0.05
results <- geneConnector(geneList=geneList,
                        networkGraph=graphReduced,
                        directed=FALSE,
                        pValueAdj="BH",
                        pValueCutoff=threshold,
                        communityMethod="ebc",
                        keepIsolatedNodes=FALSE)

# Check the p-value of the selected linker
linkerDF <- results$neighborData
linkerDF[linkerDF$pValueFDR<threshold,]

##
## The geneConnector function returns a list of data frames. 
names(results)

# plot graph with the Fruchterman-Reingold layout algorithm
plot(results$netboxCommunity,results$netboxGraph, layout=layout_with_fr) 

## ----netboxrExampleGlobalTest, eval=FALSE----------------------------------
#  ## This function will need a lot of time to complete.
#  globalTest <- globalNullModel(netboxGraph=results$netboxGraph, networkGraph=graphReduced, iterations=10, numOfGenes = 274)

## ----netboxrExampleLocalTest-----------------------------------------------
localTest <- localNullModel(netboxGraph=results$netboxGraph, iterations=1000)


## ----netboxrExampleLocalTestPlot-------------------------------------------
h<-hist(localTest$randomModularityScore,breaks=35,plot=FALSE)
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE,ylim=c(0,0.1),xlim=c(0.1,0.6), col="lightblue")
abline(v=localTest$modularityScoreObs,col="red")

## --------------------------------------------------------------------------
DT::datatable(results$moduleMembership, rownames = FALSE)

## ----netboxrEampleOutput, eval=FALSE---------------------------------------
#  # Write results for further visilaztion in the cytoscape software.
#  #
#  # network.sif file is the NetBox algorithm output in SIF format.
#  write.table(results$netboxOutput, file="network.sif", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
#  #
#  # netighborList.txt file contains the information of all neighbor nodes.
#  write.table(results$neighborData, file="neighborList.txt", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#  #
#  # community.membership.txt file indicates the identified pathway module numbers.
#  write.table(results$moduleMembership, file="community.membership.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)
#  #
#  # nodeType.txt file indicates the node is "linker" node or "candidate" node.
#  write.table(results$nodeType,file="nodeType.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

## --------------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)

module <- 6
selectedModule <- results$moduleMembership[results$moduleMembership$membership == module,]
geneList <-selectedModule$geneSymbol

# Check available ID types in for the org.Hs.eg.db annotation package
keytypes(org.Hs.eg.db)

ids <- bitr(geneList, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
head(ids)

ego <- enrichGO(gene = ids$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable = TRUE)

## --------------------------------------------------------------------------
head(ego)

## ---- fig.width=10, fig.height=5-------------------------------------------
dotplot(ego)

## ----paxtoolsr, fig.width=15, fig.height=15--------------------------------
##library(paxtoolsr)

##filename <- "PathwayCommons.8.reactome.EXTENDED_BINARY_SIF.hgnc.txt.gz"
##sif <- downloadPc2(filename, version="8")

# NOTE: Run without filename to see a list of available files
#sif <- downloadPc2()

# Filter interactions for specific types
##interactionTypes <- getSifInteractionCategories()
# filteredSif <- filterSif(sif$edges, interactionTypes=interactionTypes[["BetweenProteins"]])
##filteredSif <- filteredSif[(filteredSif$INTERACTION_TYPE %in% "in-complex-with"), ]

# Re-run NetBox algorithm with new network
##graphReduced <- networkSimplify(filteredSif, directed=FALSE)      
##geneList <- as.character(netbox2010$geneList) 

##threshold <- 0.05
##pcResults <- geneConnector(geneList=geneList,
##                           networkGraph=graphReduced,
##                           directed=FALSE,
##                           pValueAdj="BH",
##                           pValueCutoff=threshold,
##                           communityMethod="lec",
##                           keepIsolatedNodes=FALSE)

# Check the p-value of the selected linker
##linkerDF <- results$neighborData
##linkerDF[linkerDF$pValueFDR<threshold,]

# The geneConnector function returns a list of data frames. 
##names(results)

# plot graph with the Fruchterman-Reingold layout algorithm
##plot(results$netboxCommunity,results$netboxGraph, layout=layout_with_fr) 

## ----sessionInfo-----------------------------------------------------------
sessionInfo()

