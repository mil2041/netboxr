#' Union multiple data frames for making oncoprint
#' 
#' @param mat List of multiple data frames
#' @param verbose Default is TRUE
#' 
#' @return a single data frame that is union of multiple data frames
#'   
#' @examples
#' date<-'2015_08_21'
#' cancerType<-'KIRC'
#' worDir<-getwd()
#' #mutationMAF<-getMutationMAF(date,cancerType,workDir)
#' 
#' #alterationMatSingle<-unionMultipleFrames(alterationMatMultiple)
#' 
#' @concept netboxr
#' @export 
unionMultipleFrames <- function(mat, verbose = TRUE) {
    
    mat1 <- mat[[1]]
    mat2 <- mat[[2]]
    mat3 <- mat[[3]]
    
    tmpMat <- {
    }
    
    if (!is.null(nrow(mat1)) & !is.null(nrow(mat2))) {
        tmpMat <- unionTwoFrames(mat1, mat2)
    }
    
    if (is.null(nrow(mat1)) & !is.null(nrow(mat2))) {
        tmpMat <- mat2
    }
    
    if (!is.null(nrow(mat2)) & is.null(nrow(mat2))) {
        tmpMat <- mat1
    }
    
    if (is.null(nrow(mat1)) & is.null(nrow(mat2))) {
        if (is.null(nrow(mat3))) {
            alterationMat <- {
            }
        } else {
            alterationMat <- mat3
        }
    }
    
    if (!is.null(nrow(mat3))) {
        alterationMat <- unionTwoFrames(tmpMat, mat3)
    } else {
        alterationMat <- tmpMat
    }
    
    return(alterationMat)
    
}

