#' Retrieve TCGA analysis run date
#' 
#' @param latest TRUE or FALSE, if this parameter is FALSE, return all analysis date.
#' 
#' @return a vector of ananlysis date
#'   
#' @examples
#' getRunDates(latest=FALSE)
#' 
#' @concept netboxr
#' @export
#' @import XML
getRunDates <- function(latest = FALSE) {
    
    url <- "http://gdac.broadinstitute.org/runs"
    doc <- htmlTreeParse(url, useInternalNodes = TRUE)
    
    # the date between stddate and analyses may not always the same keyWord = paste('','stddata_',sep='')
    
    keyWord = paste("", "analyses_", sep = "")
    keyWord = paste("//a[contains(@href, '", keyWord, "')]", sep = "")
    plinks1 = xpathSApply(doc, keyWord, xmlValue)
    runDate <- gsub("/", "", substring(plinks1, 11))
    remove <- c("latest")
    runDate <- runDate[!runDate %in% remove]
    runDate <- runDate[rev(order(as.Date(runDate, format = "%Y_%m_%d")))]
    if (!latest) {
        return(runDate)
    } else {
        return(runDate[1])
    }
}
