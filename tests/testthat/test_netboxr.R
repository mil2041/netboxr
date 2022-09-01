# All tests are done on files in package using system.file()

test_that("sanity_check", {
  expect_equal(2 * 2, 4)
})

test_that("annotateGraph", {
  #Remove later 2
  data(netbox2010)
  interaction_type_color <- read.csv(system.file("interaction_type.color.txt", package = "netboxr"),
                                     header=TRUE, sep="\t", stringsAsFactors=FALSE)
  sifNetwork<-pathway_commons_v8_reactome$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE)
  
  geneList <- netbox2010$geneList
  results <- geneConnector(geneList = geneList, networkGraph = graphReduced, 
                           directed = FALSE, pValueAdj = "BH", pValueCutoff = 0.05, 
                           communityMethod = "lec", keepIsolatedNodes = FALSE)
  
  netboxGraphAnnotated <- annotateGraph(netboxResults = results,
                                        edgeColors = interaction_type_color,
                                        directed = TRUE,
                                        linker = TRUE)
  
  tmp_file <- read.table(system.file(file.path("test_output", "graphAnnotated.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  expect_equal(as_data_frame(netboxGraphAnnotated), tmp_file)
  })

test_that("geneConnector", {
  
  tmp_file_1 <- read.table(system.file(file.path("test_output", "network.sif"), package = "netboxr"), 
                         header=FALSE, sep="\t", stringsAsFactors=FALSE)
  tmp_file_2 <- read.table(system.file(file.path("test_output", "neighborList.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  tmp_file_3 <- read.table(system.file(file.path("test_output", "community.membership.txt"), package = "netboxr"), 
                         header=FALSE, sep="\t", stringsAsFactors=FALSE)
  tmp_file_4 <- read.table(system.file(file.path("test_output", "nodeType.txt"), package = "netboxr"), 
                         header=FALSE, sep="\t", stringsAsFactors=FALSE)
  
  data(netbox2010)
  sifNetwork<-netbox2010$network
  graphReduced <- networkSimplify(sifNetwork,directed = FALSE) 
  geneList<-as.character(netbox2010$geneList)
  
  results<-geneConnector(geneList=geneList,networkGraph=graphReduced,
                         pValueAdj='BH',pValueCutoff=0.05,
                         communityMethod='lec',keepIsolatedNodes=FALSE)
  
  expect_equal(as_data_frame(results$netboxGraph)[,1], tmp_file_1[,1])
  expect_equal(as_data_frame(results$netboxGraph)[,2], tmp_file_1[,3])
  expect_equal(as_data_frame(results$netboxGraph)[,3], tmp_file_1[,2])
  expect_equal(results$neighborData$globalDegree, tmp_file_2$globalDegree)
  expect_equal(results$moduleMembership[,1], tmp_file_3[,1])
  expect_equal(results$nodeType[,1], tmp_file_4[,1])
  expect_equal(results$nodeType[,2], tmp_file_4[,2])
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
                                iterations=5, numOfGenes = 274)
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
})

test_that("networkSimplify", {
  
  data(netbox2010)
  
  sifNetwork <- netbox2010$network
  graphReduced <- networkSimplify(sifNetwork, directed = FALSE)
  
  tmp_file <- read.table(system.file(file.path("test_output", "graphReduced.txt"), package = "netboxr"), 
                         header=TRUE, sep="\t", stringsAsFactors=FALSE)
  expect_equal(as_data_frame(graphReduced), tmp_file)
})