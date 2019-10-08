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
localNullModel <- function(netboxGraph, iterations = 30) {
  community <- edge.betweenness.community(netboxGraph)
  moduleMembership <- membership(community)
  modularityScoreObs <- modularity(netboxGraph, moduleMembership)
  # cat(sprintf('Observed Network Modularity is:%s\n',modularityScoreObs))

  modularityScore <- {
  }

  # pb <- txtProgressBar(1, iterations, style = 3)
  timeStart <- Sys.time()

  for (iter in 1:iterations) {
    # cat(sprintf('local null model iteration: %s / %s\n',iter,iterations))
    graphTmp <- rewire(netboxGraph, with = keeping_degseq(niter = vcount(netboxGraph) * 10))
    community <- edge.betweenness.community(graphTmp)
    moduleMembership <- membership(community)
    modularityScore[iter] <- modularity(graphTmp, moduleMembership)
    # cat(sprintf('Random Network Modularity is:%s\n',modularityScore[iter]))

    # if(iter %% 50 == 0) { setTxtProgressBar(pb, iter) }
  }

  Sys.time() - timeStart

  randomMean <- mean(modularityScore)
  randomSD <- sd(modularityScore)

  zScore <- (modularityScoreObs - randomMean) / randomSD
  pValueObs <- (1 - pnorm(zScore))

  cat(sprintf("###########\n"))
  cat(sprintf("Based on %s random trails\n", iterations))
  cat(sprintf("Random networks: mean modularity = %s\n", round(randomMean, digits = 3)))
  cat(sprintf("Random networks: sd modularity = %s\n", round(randomSD, digits = 3)))
  cat(sprintf("Observed network modularity is: %s\n", round(modularityScoreObs, digits = 3)))
  cat(sprintf("Observed network modularity z-score is: %s\n", round(zScore, digits = 3)))
  cat(sprintf("One-tail p-value is: %s\n", signif(pValueObs, digits = 4)))

  # d<-density(modularityScore) plot(d,xlim=c(0,1),main='Kernal Density of Local Modularity Score')
  # polygon(d, col='light blue', border='black') abline(v=modularityScoreObs,col='red')

  # h<-hist(modularityScore,breaks=35,xlim=c(0.1,0.6)) h$density = h$counts/sum(h$counts)
  # plot(h,freq=FALSE,ylim=c(0,0.1),xlim=c(0.1,0.6), col='lightblue') abline(v=modularityScoreObs,col='red')
  output <- list(
    randomModularityScore = modularityScore, randomMean = randomMean, randomSD = randomSD,
    modularityScoreObs = modularityScoreObs, zScore = zScore, pValueObs = pValueObs
  )
  return(output)
}
