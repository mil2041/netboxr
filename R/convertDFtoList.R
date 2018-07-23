#' Convert a data.frame to a list
#' 
#' @param df a data.frame
#' @param listItemNameCol the name of the column that contains the categories to be split on
#' @param listItemCol the name of the column that contains the entries to be added to list entries
#' 
#' @return a list for the given data.frame
#' 
#' @examples 
#' #df <- results$moduleMembership
#' #convertDFtoList(df, "membership", "geneSymbol")
#' 
#' @concept netboxr
#' @export
convertDFtoList <- function(df, listItemNameCol, listItemCol) {
  x <- list()
  
  for(i in 1:max(df[, listItemNameCol])) {
    tmp <- as.vector(df[df[, listItemNameCol] == i, listItemCol])
    
    x[[as.character(i)]] <- tmp
  }
  
  return(x)
}
