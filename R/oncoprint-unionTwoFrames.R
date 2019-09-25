#' Union data frames for making oncoprint
#' 
#' @param mat1 data frame 1
#' @param mat2 data frame 2
#' @param verbose Default is TRUE
#' 
#' @return a union data frame of mat1 and mat2
#'   
#' @examples
#' date<-'2015_08_21'
#' cancerType<-'KIRC'
#' worDir<-getwd()
#' #mutationMAF<-getMutationMAF(date,cancerType,workDir)
#' 
#' #alterationMat<-unionTwoFrames(memoMutationMat,memoCNCallMat)
#' 
#' @concept netboxr
#' @export 
#' @importFrom reshape2 acast
unionTwoFrames <- function(mat1, mat2, verbose = TRUE) {
    
    mat1 <- as.matrix(mat1)
    mat2 <- as.matrix(mat2)
    mat1[is.na(mat1)] <- 0
    mat2[is.na(mat2)] <- 0
    
    mat1ColName <- substr(colnames(mat1), 1, 12)
    mat2ColName <- substr(colnames(mat2), 1, 12)
    
    unionColName <- union(mat1ColName, mat2ColName)
    unionRowName <- union(rownames(mat1), rownames(mat2))
    
    # initialDF<-{}
    
    # for(i in 1:length(unionRowName)){ a<-rep(unionRowName[i],length(unionColName)) b<-unionColName
    # c<-rep(0,length(unionColName)) tmp1<-data.frame(a,b,c) tmp1$index<-paste(tmp1$a,tmp1$b,sep='@')
    # initialDF<-rbind(initialDF,tmp1) }
    
    mm <- matrix(0, length(unionRowName), length(unionColName))
    colnames(mm) <- unionColName
    rownames(mm) <- unionRowName
    tt <- cbind(expand.grid(dimnames(mm)), value = as.vector(mm))
    tt$index <- paste(tt[, 1], tt[, 2], sep = "@")
    initialDF <- tt
    
    colnames(initialDF) <- c("rowName", "colName", "value", "index")
    
    mat1DF <- cbind(expand.grid(dimnames(mat1)), value = as.vector(mat1))
    colnames(mat1DF) <- c("rowName", "colName", "value")
    mat1DF$index <- paste(mat1DF[, 1], substr(mat1DF[, 2], 1, 12), sep = "@")
    
    mat2DF <- cbind(expand.grid(dimnames(mat2)), value = as.vector(mat2))
    colnames(mat2DF) <- c("rowName", "colName", "value")
    mat2DF$index <- paste(mat2DF[, 1], substr(mat2DF[, 2], 1, 12), sep = "@")
    
    tmp1 <- merge(initialDF, mat1DF, by = c("index"), all = TRUE)
    tmp2 <- merge(tmp1, mat2DF, by = c("index"), all = TRUE)
    tmp3 <- data.frame(tmp2$index, tmp2[7], tmp2[10])
    tmp3[is.na(tmp3)] <- 0
    
    tmp3$catStr <- rep(0, nrow(tmp3))
    tmp3$catStr <- paste(tmp3[, 2], tmp3[, 3], sep = "")
    tmp3$catStr <- gsub("0", "", tmp3$catStr)
    tmp3[tmp3$catStr == "", ]$catStr <- 0
    colnames(tmp3) <- c("index", "v1", "v2", "catStr")
    
    tmp4 <- data.frame(do.call("rbind", strsplit(as.character(tmp3$index), "@", fixed = TRUE)))
    tmp5 <- data.frame(tmp4, tmp3$catStr)
    colnames(tmp5) <- c("rowName", "colName", "catStr")
    tmp5 <- data.frame(lapply(tmp5, as.character), stringsAsFactors = FALSE)
    alterationMat <- acast(tmp5, rowName ~ colName, value.var = "catStr")
    alterationMat[alterationMat == 0] <- NA
    
    alterationMat <- as.data.frame(alterationMat, stringsAsFactors = FALSE)
    
    return(alterationMat)
    
}

