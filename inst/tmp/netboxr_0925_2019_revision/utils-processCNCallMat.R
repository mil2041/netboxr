#' Prcoess GISTIC2 copy number call from TCGA study
#' 
#' @param CNCallMat GISTIC2 copy number call data frame
#' @param verbose output detailed information [ TRUE, FALSE ] 
#' 
#' @return a data frame of copy number call from GISTIC2 reuslt
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-'KIRC'
#' #worDir<-getwd()
#' 
#' #gistic2CNCallMat<-processCNCallMat(gistic2CNCallMat)
#' 
#' @concept netboxr
#' @export
processCNCallMat <- function(CNCallMat, verbose = TRUE) {
    
    cat(sprintf("Pre-processing GISTIC2 CN call matrix\n"))
    
    # entrezID<-CNCallMat[,2] geneSymbol<-CNCallMat[,1]
    # CNCallgeneSymbolTable<-data.frame(geneSymbol,entrezID) CNCallMat<-CNCallMat[,4:ncol(CNCallMat)]
    # rownames(CNCallMat)<-entrezID sampleIDs<-gsub('-','.',colnames(CNCallMat))
    
    geneSymbol <- CNCallMat[, 1]
    locusID <- CNCallMat[, 2]
    CNCallgeneSymbolTable <- data.frame(geneSymbol, locusID, stringsAsFactors = FALSE)
    CNCallMat <- CNCallMat[, 4:ncol(CNCallMat)]
    rownames(CNCallMat) <- geneSymbol
    sampleIDs <- colnames(CNCallMat)
    
    ### code from RTCGAToolbox.R sampleIDs<-as.character(id_df[,1])
    samplesDat <- data.frame(matrix(nrow = length(sampleIDs), ncol = 7))
    rownames(samplesDat) <- sampleIDs
    for (j in 1:length(sampleIDs)) {
        tmpRow <- unlist(strsplit(sampleIDs[j], split = "-", fixed = TRUE))
        samplesDat[sampleIDs[j], ] <- tmpRow
    }
    sampleIDs1 <- as.character(samplesDat[, 4])
    sampleIDs1 <- substr(sampleIDs1, 1, nchar(sampleIDs1) - 1)
    sampleIDs1 <- as.numeric(sampleIDs1)
    normalSamples <- rownames(samplesDat)[sampleIDs1 < 20 & sampleIDs1 > 9]
    tumorSamples <- rownames(samplesDat)[sampleIDs1 < 2]
    
    if (length(tumorSamples) == 0) {
        tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
    }
    
    #### CNCallSamples<-tumorSamples
    
    if (length(normalSamples) > 0) {
        # CNCallSamplesNormal<-normalSamples CNCallMatNormal<-CNCallMat[,gsub('\\.','-',normalSamples)]
        CNCallMatNormal <- CNCallMat[, normalSamples]
    } else {
        CNCallMatNormal <- NULL
    }
    
    # CNCallMatTumor<-CNCallMat[,gsub('\\.','-',tumorSamples)]
    CNCallMatTumor <- CNCallMat[, tumorSamples]
    
    cat(sprintf("The file contains %s tumor samples and %s normal samples\n", length(tumorSamples), length(normalSamples)))
    
    CNMat <- list(tumorMat = CNCallMatTumor, normalMat = CNCallMatNormal, geneSymbolTable = CNCallgeneSymbolTable)
    
    return(CNMat)
    
}
