#' Prcoess GISTIC2 relative copy number value from TCGA study
#' 
#' @param relativeCNinGISTIC GISTIC2 relative copy number value data frame
#' @param verbose output detailed information [ TRUE, FALSE ] 
#' 
#' @return a data frame of absolute copy number value from GISTIC2 reuslt
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-'KIRC'
#' #worDir<-getwd()
#' 
#' #absoluteCNValueMat<-processCNValueMat(GISTIC2CNValueMat)
#' 
#' @concept netboxr
#' @export
processCNValueMat <- function(relativeCNinGISTIC, verbose = TRUE) {
    
    cat(sprintf("Pre-processing GISTIC2 CN value matrix\n"))
    # setwd('~/work/Sander_lab/TCGA_data/BRCA/TCGA/CopyNumber/GISTIC2') CN value in unit of ( CN - 2 )
    # relativeCNinGISTIC<-read.table('all_data_by_genes.txt.gz',sep='\t',header=TRUE,quote=,'',comment='')
    
    locusID <- relativeCNinGISTIC[, 2]
    geneSymbol <- relativeCNinGISTIC[, 1]
    CNValuegeneSymbolTable <- data.frame(geneSymbol, locusID, stringsAsFactors = FALSE)
    
    absoluteCNMat <- relativeCNinGISTIC[, 4:ncol(relativeCNinGISTIC)] + 2
    # CNValueMat<-data.frame(absoluteCNMat)
    CNValueMat <- absoluteCNMat
    rownames(CNValueMat) <- geneSymbol
    
    sampleIDs <- colnames(CNValueMat)
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
    
    ## the TCGA guide said tumor samples are ID<10, but we only use 01A
    tumorSamples <- rownames(samplesDat)[sampleIDs1 < 2]
    
    if (length(tumorSamples) == 0) {
        tumorSamples <- rownames(samplesDat)[sampleIDs1 < 10]
    }
    
    #### CNValueSamplesTumor<-tumorSamples
    
    if (length(normalSamples) > 0) {
        # CNValueSamplesNormal<-normalSamples
        CNValueMatNormal <- CNValueMat[, normalSamples]
    } else {
        CNValueMatNormal <- NULL
    }
    
    CNValueMatTumor <- CNValueMat[, tumorSamples]
    
    cat(sprintf("The file contains %s tumor samples and %s normal samples\n", length(tumorSamples), length(normalSamples)))
    
    CNMat <- list(tumorMat = CNValueMatTumor, normalMat = CNValueMatNormal, geneSymbolTable = CNValuegeneSymbolTable)
    
    return(CNMat)
    
}
