#' Make binary mutation matrix for oncoprint
#' 
#' @param mutationMAF mutation MAF data frame
#' @param selectedSampleId subset of sample id
#' @param geneList Where to save the processed result
#' @param verbose 
#' 
#' @return a binary mutation data frame
#'   
#' @examples
#' date<-'2015_08_21'
#' cancerType<-'KIRC'
#' worDir<-getwd()
#' #mutationMAF<-getMutationMAF(date,cancerType,workDir)
#' 
#' #mutationMat<-makeMutationMatrix(mutationMAF,selectedSampleId,geneList)
#' 
#' @concept netboxr
#' @export 
#' @importFrom reshape2 acast
makeMutationMatrix <- function(mutationMAF, selectedSampleId, geneList, verbose) {
    
    tmpMAF <- mutationMAF[which(substr(mutationMAF$patientID, 1, 12) %in% substr(selectedSampleId, 1, 12)), 
        ]
    overlapSampleID <- intersect(unique(substr(mutationMAF$patientID, 1, 12)), unique(substr(selectedSampleId, 
        1, 12)))
    cat(sprintf("%s samples overlap with %s samples from query\n", length(overlapSampleID), length(selectedSampleId)))
    
    overlapGene <- intersect(unique(tmpMAF$Hugo_Symbol), geneList)
    cat(sprintf("%s genes overlap with %s genes from query\n", length(overlapGene), length(geneList)))
    
    keep <- tmpMAF$Hugo_Symbol %in% geneList
    
    if (!is.null(keep)) {
        tmpMAF <- tmpMAF[keep, ]
    } else {
        tmpMAF <- {
        }
    }
    
    if (nrow(tmpMAF) > 0) {
        
        tmpDF <- data.frame(tmpMAF$Hugo_Symbol, tmpMAF$patientID, stringsAsFactors = FALSE)
        colnames(tmpDF) <- c("Hugo_Symbol", "patientID")
        tmpDF$check <- paste(tmpDF[, 1], tmpDF[, 2], sep = "@")
        tmpDF <- tmpDF[order(tmpDF$check), ]
        removed <- duplicated(tmpDF$check)
        tmpDF <- tmpDF[!removed, ]
        # c<-c[,-3] c$value<-rep('1',nrow(c))
        
        mutationGeneInMemo <- unique(tmpMAF$Hugo_Symbol)
        keep <- substr(mutationMAF$patientID, 1, 12) %in% substr(selectedSampleId, 1, 12)
        mutationSamplesID <- unique(mutationMAF$patientID[keep])
        mutationMat <- matrix(0, ncol = length(mutationSamplesID), nrow = length(mutationGeneInMemo))
        colnames(mutationMat) <- mutationSamplesID
        rownames(mutationMat) <- mutationGeneInMemo
        
        initialMat <- cbind(expand.grid(dimnames(mutationMat)), value = as.vector(mutationMat))
        colnames(initialMat) <- c("Hugo_Symbol", "patientID", "value")
        initialMat$check <- paste(initialMat[, 1], initialMat[, 2], sep = "@")
        
        replaceRow <- initialMat$check %in% tmpDF$check
        initialMat[replaceRow, ]$value <- 1
        initialMat <- initialMat[, -4]
        
        # use reshape2 package to make mutation matrix
        mutationMat <- acast(initialMat, Hugo_Symbol ~ patientID, value.var = "value")
        mutationMat <- as.data.frame(mutationMat)
    } else {
        mutationMat <- {
        }
    }
    
    memoMutationMat <- mutationMat
    
    memoMutationMat[memoMutationMat == (0)] <- NA
    memoMutationMat[memoMutationMat == (1)] <- "MUT;"
    
    return(memoMutationMat)
    
}
