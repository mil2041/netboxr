# All tests are done on files in package using system.file()

test_that("sanity_check", {
  expect_equal(2 * 2, 4)
})

test_that("annotateGraph", {
  data(netbox2010)
  interaction_type_color <- read.csv(system.file("interaction_type.color.txt", package = "netboxr"),
                                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE)
  
  geneList <- netbox2010$geneList
  results <- geneConnector(geneList = geneList, networkGraph = graphReduced, 
                           directed = FALSE, pValueAdj = "BH", pValueCutoff = 0.05, 
                           communityMethod = "lec", keepIsolatedNodes = FALSE)
  
  library(RColorBrewer)
  edges <- results$netboxOutput
  interactionType<-unique(edges[,2])
  interactionTypeColor<-brewer.pal(length(interactionType),name="Spectral")
  
  edgeColors<-data.frame(interactionType,interactionTypeColor,stringsAsFactors = FALSE)
  colnames(edgeColors)<-c("INTERACTION_TYPE","COLOR")
  
  netboxGraphAnnotated <- annotateGraph(netboxResults = results,
                                        edgeColors = edgeColors,
                                        directed = TRUE,
                                        linker = TRUE)
  
  tmp_file <- read.table(system.file(file.path("test_output", "graphAnnotated.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char = "")
  expect_equal(as_data_frame(netboxGraphAnnotated), tmp_file)
})

test_that("geneConnector", {
  
  tmp_file_1 <- read.table(system.file(file.path("test_output", "network.sif"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  tmp_file_2 <- read.table(system.file(file.path("test_output", "neighborList.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  tmp_file_3 <- read.table(system.file(file.path("test_output", "community.membership.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  tmp_file_4 <- read.table(system.file(file.path("test_output", "nodeType.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='ebc',keepIsolatedNodes=FALSE)
  
  expect_equal(results$netboxOutput, tmp_file_1)
  expect_equal(results$neighborData$globalDegree, tmp_file_2$globalDegree)
  expect_equal(results$moduleMembership$geneSymbol, tmp_file_3$geneSymbol)
  expect_equal(results$nodeType, tmp_file_4)
})

test_that("globalNullModel", {
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='ebc',keepIsolatedNodes=FALSE)
  
  globalTest <- globalNullModel(netboxGraph=results$netboxGraph, 
                                networkGraph=graphReduced, 
                                iterations=10, numOfGenes = 274)
  
  expect_equal(length(globalTest), 3)
  expect_gt(globalTest$globalNodesResult$pValueNodes, 0.09)
  expect_gt(globalTest$globalEdgesResult$pValueEdges, 0.09)
})

test_that("localNullModel", {
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='ebc',keepIsolatedNodes=FALSE)
  
  localTest <- localNullModel(netboxGraph=results$netboxGraph, iterations=10)
  
  load(system.file(file.path("test_output", "localTest_test.rda"), package = "netboxr"))
  expect_equal(localTest_test$modularityScoreObs, localTest$modularityScoreObs)
})

test_that("networkSimplify", {
  
  data(netbox2010)
  
  sifNetwork <- netbox2010$network
  graphReduced <- networkSimplify(sifNetwork, directed = FALSE)
  
  tmp_file <- read.table(system.file(file.path("test_output", "graphReduced.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  expect_equal(as_data_frame(graphReduced), tmp_file)
})
