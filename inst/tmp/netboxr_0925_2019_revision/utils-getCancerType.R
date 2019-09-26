#' Retrieve TCGA analysis summary
#' 
#' @param date 
#' @param cancerType 
#' 
#' @return a data frame of cancer study summary
#'   
#' @examples
#' date<-'2015_08_21'
#' getCancerType(date)
#' 
#' @concept netboxr
#' @export
#' @import XML
getCancerType <- function(date = "last", cancerType = "all") {
    
    
    if (date %in% "last") {
        date <- getRunDates(last = TRUE)
    }
    
    cat(sprintf("Retreive data from %s\n", date))
    
    url <- "http://gdac.broadinstitute.org/runs"
    url <- paste(url, "/stddata__", date, sep = "")
    url <- paste(url, "/samples_report", sep = "")
    doc <- htmlTreeParse(url, useInternalNodes = TRUE)
    
    keyWord = paste("", "sample_counts", sep = "")
    keyWord = paste("//a[contains(@href, '", keyWord, "')]", sep = "")
    plinks1 = xpathSApply(doc, keyWord, xmlAttrs)
    data_fn <- plinks1[2]
    
    url <- "http://gdac.broadinstitute.org/runs"
    url <- paste(url, "/stddata__", date, sep = "")
    url <- paste(url, "/samples_report", sep = "")
    url <- paste(url, "/", data_fn, sep = "")
    
    data <- read.table(url, header = TRUE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE, 
        stringsAsFactors = FALSE)
    
    data <- data[!data$Cohort %in% "Totals", ]
    
    if (cancerType %in% "all") {
        detailedSummaryLine <- grep("-", data$Cohort)
        keep <- setdiff(1:nrow(data), detailedSummaryLine)
        return(data[keep, ])
    }
    
    if (cancerType %in% "NT") {
        keep <- grep(cancerType, data$Cohort)
        return(data[keep, ])
    }
    
    if (cancerType %in% "TP") {
        keep <- grep(cancerType, data$Cohort)
        return(data[keep, ])
    }
    
    if (cancerType %in% "full") {
        return(data)
    }
    
}





