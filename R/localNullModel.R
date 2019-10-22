#' @title Generate local null model p-value
#'
#' @description
#' This function keeps the number of connections of each nodes in the graph but
#' it rewires the partners of connections and produces  modularity score. When
#' it repeats multiple time, a modularity score distribution will be used to
#' produce netbox loacl p-value.
#'
#' @param netboxGraph A vector containing candidate gene list
#'
#' @param iterations TRUE of FALSE
#'
#' @return a list with four lists (i.e. netboxOutput, nodeType, moduleMembership, neighborData)
#'
#' @author Eric Minwei Liu, \email{emliu.research@gmail.com}
#'
#' @examples
#' data(netbox2010)
#'
#' # geneList<-netbox2010$geneList
#' # sifNetwork<-netbox2010$network
#' # result<-geneConnector(geneList=geneList,sifNetwork=network,pValueAdj='BH',
#' #           pValueCutoff=0.05,communityMethod='lec',keepIsolatedNodes=FALSE)
#'
#' # names(result)
#' @concept netboxr
#' @export
#' @import igraph
#' @importFrom stats sd
#' @importFrom stats pnorm
localNullModel <- function(netboxGraph, iterations = 30) {
  community <- edge.betweenness.community(netboxGraph)
  moduleMembership <- membership(community)
  modularityScoreObs <- modularity(netboxGraph, moduleMembership)

  modularityScore <- {
  }

  timeStart <- Sys.time()

  for (iter in seq_len(iterations)) {
    # message(sprintf('local null model iteration: %s / %s\n',iter,iterations))
    graphTmp <- rewire(netboxGraph, with = keeping_degseq(niter = vcount(netboxGraph) * 10))
    community <- edge.betweenness.community(graphTmp)
    moduleMembership <- membership(community)
    modularityScore[iter] <- modularity(graphTmp, moduleMembership)
    
  }

  Sys.time() - timeStart

  randomMean <- mean(modularityScore)
  randomSD <- sd(modularityScore)

  zScore <- (modularityScoreObs - randomMean) / randomSD
  pValueObs <- (1 - pnorm(zScore))

  message(sprintf("###########\n"))
  message(sprintf("Based on %s random trails\n", iterations))
  message(sprintf("Random networks: mean modularity = %s\n", round(randomMean, digits = 3)))
  message(sprintf("Random networks: sd modularity = %s\n", round(randomSD, digits = 3)))
  message(sprintf("Observed network modularity is: %s\n", round(modularityScoreObs, digits = 3)))
  message(sprintf("Observed network modularity z-score is: %s\n", round(zScore, digits = 3)))
  message(sprintf("One-tail p-value is: %s\n", signif(pValueObs, digits = 4)))

  output <- list(
    randomModularityScore = modularityScore, randomMean = randomMean, randomSD = randomSD,
    modularityScoreObs = modularityScoreObs, zScore = zScore, pValueObs = pValueObs
  )
  return(output)
}
