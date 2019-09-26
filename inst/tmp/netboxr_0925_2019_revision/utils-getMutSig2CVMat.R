#' Retrieve significant mutation list from TCGA study
#' 
#' @param date Default is 'latest' 
#' @param cancerType Abbreviation in the TCGA study ['KIRC','GBM']
#' @param selectedSampleId Subset of TCGA sample id
#' @param workDir Where to save the processed result
#' 
#' @return a data frame of MutSig2CV result
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-'KIRC'
#' #worDir<-getwd()
#' 
#' #mutSig2CVMat<-getMutSig2CVMat(date,cancerType,workDir)
#' 
#' @concept netboxr
#' @export 
#' @import XML
getMutSig2CVMat <- function(date = "last", cancerType = "BRCA", workDir, verbose = TRUE) {
    
    destDir <- workDir
    
    if (!file.exists(file.path(destDir, "sig_genes.txt"))) {
        
        # date<-'2014_12_06'
        if (date %in% "last") {
            date <- getRunDates(last = TRUE)
        }
        
        cat(sprintf("Retreive %s MutSig2CV data from %s\n", cancerType, date))
        
        url <- "http://gdac.broadinstitute.org/runs"
        url <- paste(url, "/analyses__", date, sep = "")
        url <- paste(url, "/data/", cancerType, "/", substr(date, 1, 4), substr(date, 6, 7), substr(date, 
            9, 10), sep = "")
        doc <- htmlTreeParse(url, useInternalNodes = TRUE)
        
        keyWord = paste("", "MutSigNozzleReport2CV.Level_4", sep = "")
        keyWord = paste("//a[contains(@href, '", keyWord, "')]", sep = "")
        plinks = xpathSApply(doc, keyWord, xmlValue)
        plinks = plinks[grepl("*.MutSigNozzleReport2CV.Level_4.*.tar[.]gz$", plinks)]
        plinks <- gsub("\\s", "", plinks)
        
        
        if (!file.exists(destDir)) {
            dir.create(file.path(destDir, cancerType), recursive = TRUE)
        }
        
        download_link = paste(url, "/", plinks, sep = "")
        download.file(url = download_link, destfile = file.path(destDir, "MutSig2CV.tar.gz"), method = "auto", 
            quiet = FALSE, mode = "w")
        fileList <- untar(file.path(destDir, "MutSig2CV.tar.gz"), list = TRUE)
        grepSearch = paste("sig_genes.txt", sep = "")
        fileList = fileList[grepl(grepSearch, fileList)]
        untar(file.path(destDir, "MutSig2CV.tar.gz"), files = fileList)
        
        file.rename(from = fileList, to = file.path(destDir, "sig_genes.txt"))
        file.remove(file.path(destDir, "MutSig2CV.tar.gz"))
        delFodler <- file.path(getwd(), strsplit(fileList, "/")[[1]][1])
        message(delFodler)
        unlink(delFodler, recursive = TRUE)
    }
    
    cat(sprintf("Loading MutSig2CV: sig_genes.txt\n"))
    
    tmpMat <- read.table(file = file.path(destDir, "sig_genes.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    return(tmpMat)
    
    
}
