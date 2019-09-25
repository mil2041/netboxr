#' Retrieve annotated mutation MAF file from TCGA study
#' 
#' @param date Default is 'latest' 
#' @param cancerType Abbreviation in the TCGA study ['KIRC','GBM']
#' @param workDir Where to save the processed result
#' @param verbose 
#' 
#' @return a data frame of annotated mutation MAF file
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-'KIRC'
#' #worDir<-getwd()
#' 
#' #mutationMAF<-getMutationMAF(date,cancerType,workDir)
#' 
#' @concept netboxr
#' @export
#' @import XML
#' @importFrom R.utils gzip
getMutationMAF <- function(date = "last", cancerType = "BRCA", workDir, verbose = TRUE) {
    
    destDir <- workDir
    
    fileName <- paste(cancerType, "-TP.final_analysis_set.maf.gz", sep = "")
    fileName <- file.path(destDir, fileName)
    if (!file.exists(fileName)) {
        
        # date<-'2014_12_06'
        if (date %in% "last") {
            date <- getRunDates(last = TRUE)
        }
        
        cat(sprintf("Retreive %s MutSig2CV MAF data from %s\n", cancerType, date))
        # http://gdac.broadinstitute.org/runs/analyses__2015_08_21/reports/cancer/PRAD-TP/MutSigNozzleReport2CV/PRAD-TP.final_analysis_set.maf
        url <- "http://gdac.broadinstitute.org/runs"
        url <- paste(url, "/analyses__", date, sep = "")
        url <- paste(url, "/reports/cancer/", cancerType, "-TP/MutSigNozzleReport2CV/", sep = "")
        fileName <- paste(cancerType, "-TP.final_analysis_set.maf", sep = "")
        url <- paste(url, fileName, sep = "")
        
        
        if (!file.exists(destDir)) {
            dir.create(file.path(workDir, cancerType), recursive = TRUE)
        }
        
        download_link = paste(url, sep = "")
        download.file(url = download_link, destfile = file.path(destDir, fileName), method = "auto", quiet = FALSE, 
            mode = "w")
        gzip(file.path(destDir, fileName))
    }
    
    cat(sprintf("Loading %s-TP.final_analysis_set.maf.gz\n", cancerType))
    fileName <- paste(cancerType, "-TP.final_analysis_set.maf.gz", sep = "")
    tmpMat <- read.table(file = file.path(destDir, fileName), header = TRUE, sep = "\t", fill = TRUE, quote = NULL, 
        comment = "", stringsAsFactors = FALSE)
    
    return(tmpMat)
    
}
