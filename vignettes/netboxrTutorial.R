## ----knitrSetup, include=FALSE-------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", fig.width=12, fig.height=12, tidy=TRUE)

## ----style, include=FALSE, echo=FALSE, results='asis'--------------------
BiocStyle::markdown()

## ----installNetBoxr, eval=FALSE------------------------------------------
#  library(remotes)
#  install_github(repo="mil2041/netboxr", ref="master", build_vignette=TRUE)

## ----loadLibrary, message=FALSE, warning=FALSE---------------------------
library(netboxr)

## ----searchHelp, eval=FALSE, tidy=FALSE----------------------------------
#  help.search("netboxr")

## ----showHelp, eval=FALSE, tidy=FALSE------------------------------------
#  help(geneConnector)
#  ?geneConnector

## ----netboxrExampleNetwork-----------------------------------------------
data(netbox2010)
sifNetwork <- netbox2010$network
graphReduced <- networkSimplify(sifNetwork,directed = FALSE)      

## ----netboxrExampleGene--------------------------------------------------
geneList <- as.character(netbox2010$geneList) 
length(geneList)

## ----netboxrExampleGeneConnector, fig.width=12, fig.height=12------------
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

## ----netboxrExampleGlobalTest, eval=FALSE--------------------------------
#  ## This function will need a lot of time to complete.
#  globalTest <- globalNullModel(netboxGraph=results$netboxGraph, networkGraph=graphReduced, iterations=10, numOfGenes = 274)

