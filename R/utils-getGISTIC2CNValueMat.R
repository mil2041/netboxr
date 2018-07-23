#' Retrieve GISTIC2 all_data_by_genes.txt file from TCGA study
#' 
#' @param date Default is "latest" 
#' @param cancerType Abbreviation in the TCGA study ["KIRC","GBM"]
#' @param workDir Where to save the processed result
#' @param verbose output detailed information [ TRUE, FALSE ] 
#' 
#' @return a data frame of relative copy number value from GISTIC2 reuslt
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' #worDir<-getwd()
#' 
#' #gistic2CNValueMat<-getGISTIC2CNValueMat(date,cancerType,workDir)
#' 
#' @concept netboxr
#' @export
#' @import XML
#' @importFrom data.table fread
getGISTIC2CNValueMat<-function(date="last",cancerType,workDir){

    #http://gdac.broadinstitute.org/runs/analyses__2014_10_17/data/BRCA-TP/20141017/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz

    destDir<-workDir  
  
    if( !file.exists(paste(destDir,"/","all_data_by_genes.txt",sep="")) ) {

    #date<-"2014_12_06"
    if( date %in% "last"){
      date<-getRunDates(last=TRUE)
    }
    
    cat(sprintf("Retreive %s GISTIC2 data from %s\n",cancerType,date))
    
    url<-"http://gdac.broadinstitute.org/runs"
    url<-paste(url,"/analyses__",date,sep="")
    url<-paste(url,"/data/",cancerType,"/",substr(date,1,4),substr(date,6,7),substr(date,9,10),sep="")
    doc<-htmlTreeParse(url,useInternalNodes=T)
    
    keyWord = paste("","CopyNumber_Gistic2.Level_4",sep="")
    keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
    plinks = xpathSApply(doc, keyWord, xmlValue)
    plinks = plinks[grepl("*.CopyNumber_Gistic2.Level_4.*.tar[.]gz$",plinks)]
    plinks<-gsub("\\s","",plinks)
    
    
    if( !file.exists(destDir) ){
      dir.create(file.path(destDir,cancerType),recursive=TRUE)
    }
    
    download_link = paste(url,"/",plinks,sep="")
    download.file(url=download_link,destfile=file.path(destDir,"GISTIC2.tar.gz"),method="auto",quiet = FALSE, mode = "w")
    fileList <- untar(file.path(destDir,"GISTIC2.tar.gz"),list=TRUE)
    grepSearch = paste("all_data_by_genes.txt",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    untar(file.path(destDir,"GISTIC2.tar.gz"),files=fileList)
    
    file.rename(from=fileList,to=file.path(destDir,"all_data_by_genes.txt"))
    file.remove(file.path(destDir,"GISTIC2.tar.gz"))
    delFodler <- file.path(getwd(),strsplit(fileList,"/")[[1]][1])
    message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    }
    
    cat(sprintf("Loading GISTIC2: all_data_by_genes.txt\n"))
    
    tmpMat<-fread(file.path(destDir,"all_data_by_genes.txt"),sep="\t",header=TRUE,stringsAsFactors=FALSE,data.table=FALSE)
    #tmpMat<-read.table(file=paste(destDir,"/","sig_genes.txt",sep=""),header=TRUE,sep="\t")
        
    return(tmpMat)
    
    
}
