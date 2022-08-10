# All tests are done on files in package using system.file()

test_that("sanity_check", {
  expect_equal(2 * 2, 4)
})

test_that("gene_connector_ebc", {
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='ebc',keepIsolatedNodes=FALSE)
  
  expect_equal(names(results), c("netboxGraph", "netboxCommunity", "netboxOutput", 
                                 "nodeType", "moduleMembership", "neighborData"))
})

test_that("global_null_model_ebc", {
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='ebc',keepIsolatedNodes=FALSE)
  
  globalTest <- globalNullModel(netboxGraph=results$netboxGraph, networkGraph=graphReduced, 
                                iterations=10, numOfGenes = 274)
  
  expect_equal(names(globalTest), c("globalNull", "globalNodesResult", "globalEdgesResult"))
})

test_that("local_null_model_ebc", {
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='ebc',keepIsolatedNodes=FALSE)
  
  localTest <- localNullModel(netboxGraph=results$netboxGraph, iterations=10)
  
  expect_equal(names(localTest), c("randomModularityScore", "randomMean", "randomSD",
                                   "modularityScoreObs", "zScore", "pValueObs"))
})

test_that("gene_connector_lec", {
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='lec',keepIsolatedNodes=FALSE)
  
  expect_equal(names(results), c("netboxGraph", "netboxCommunity", "netboxOutput", 
                                 "nodeType", "moduleMembership", "neighborData"))
})

test_that("global_null_model_lec", {
  
  globalTest <- globalNullModel(netboxGraph=results$netboxGraph, networkGraph=graphReduced, 
                                iterations=10, numOfGenes = 274)
  
  expect_equal(names(globalTest), c("globalNull", "globalNodesResult", "globalEdgesResult"))
})
