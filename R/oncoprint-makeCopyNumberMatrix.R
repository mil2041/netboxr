#' Make binary mutation matrix for oncoprint
#' 
#' @param CNCallMat GISTIC2 copy number call data frame
#' @param selectedSampleId subset of sample id
#' @param geneList Where to save the processed result
#' @param verbose 
#' 
#' @return a binary copy number call data frame
#'   
#' @examples
#' date<-'2015_08_21'
#' cancerType<-'KIRC'
#' worDir<-getwd()
#' #mutationMAF<-getMutationMAF(date,cancerType,workDir)
#' 
#' #copyNumberMat<-makeCopyNumberMatrix(CNCallMat,selectedSampleId,geneList)
#' 
#' @concept netboxr
#' @export 
#' @importFrom reshape2 acast
makeCopyNumberMatrix <- function(CNCallMat, selectedSampleId, geneList, verbose = TRUE) {
    
    tmpCNCallMat <- CNCallMat[, which(substr(colnames(CNCallMat), 1, 12) %in% substr(selectedSampleId, 
        1, 12))]
    overlapSampleID <- intersect(unique(substr(colnames(tmpCNCallMat), 1, 12)), unique(substr(selectedSampleId, 
        1, 12)))
    cat(sprintf("%s samples overlap with %s samples from query\n", length(overlapSampleID), length(selectedSampleId)))
    
    overlapGene <- intersect(unique(rownames(tmpCNCallMat)), geneList)
    cat(sprintf("%s genes overlap with %s genes from query\n", length(overlapGene), length(geneList)))
    
    keep <- rownames(tmpCNCallMat) %in% geneList
    
    if (!is.null(keep)) {
        tmpCNCallMat <- tmpCNCallMat[keep, ]
    } else {
        tmpCNCallMat <- {
        }
    }
    
    memoCNCallMat <- tmpCNCallMat
    
    memoCNCallMat[memoCNCallMat == (1)] <- 0
    memoCNCallMat[memoCNCallMat == (-1)] <- 0
    memoCNCallMat[memoCNCallMat == (0)] <- NA
    memoCNCallMat[memoCNCallMat == (2)] <- "AMP;"
    memoCNCallMat[memoCNCallMat == (-2)] <- "HOMDEL;"
    
    return(memoCNCallMat)
    
}
