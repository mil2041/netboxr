#' Retrieve RNASeqV2 raw counts from TCGA study
#' 
#' @param date Default is "latest" 
#' @param cancerType Abbreviation in the TCGA study ["KIRC","GBM"]
#' @param workDir Where to save the processed result
#' @param verbose output detailed information [ TRUE, FALSE ] 
#' 
#' @return a data frame of copy number call from GISTIC2 reuslt
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerType<-"KIRC"
#' #worDir<-getwd()
#' 
#' #rawCountsMat<-getRNASeqV2RawCounts(date=date,cancerType=cancerType,workDir=filePath)
#' 
#' @concept netboxr
#' @export
#' @import XML
#' @importFrom data.table fread
getRNASeqV2RawCounts<-function(date="last",cancerType,workDir,verbose=TRUE){

    #http://gdac.broadinstitute.org/runs/analyses__2014_10_17/data/BRCA-TP/20141017/gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2014101700.0.0.tar.gz

    destDir<-workDir
    tmpMat<-{}
    
    if( !file.exists(file.path(destDir,"RSEM_gene_combined.Rd")) ) {
    
    #date<-"2014_12_06"
    if( date %in% "last"){
      date<-getRunDates(last=TRUE)
    }
    
    cat(sprintf("Retreive %s RSEM raw counts data from %s\n",cancerType,date))
    
    url<-"http://gdac.broadinstitute.org/runs"
    url<-paste(url,"/stddata__",date,sep="")
    url<-paste(url,"/data/",cancerType,"/",substr(date,1,4),substr(date,6,7),substr(date,9,10),sep="")
    doc<-htmlTreeParse(url,useInternalNodes=T)
    
    keyWord = paste("","RSEM_genes__data.Level_3",sep="")
    keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
    plinks = xpathSApply(doc, keyWord, xmlValue)
    plinks = plinks[grepl("*.*RSEM_genes__data.Level_3.*.tar[.]gz$",plinks)]
    plinks<-gsub("\\s","",plinks)
    #remove FFPE sample match
    tmp<-grep("FFPE",plinks)
    if(length(tmp)>0){
      plinks<-plinks[-tmp]
    }
      
    if( !file.exists(paste(destDir,sep="/")) ){
      dir.create(paste(workDir,"/",cancerType,sep=""),recursive=TRUE)
    }
    
    for (i in 1:length(plinks)){
     
    download_link = paste(url,plinks[i],sep="/")
    download.file(url=download_link,destfile=file.path(destDir,"RSEM_genes_data.tar.gz"),method="auto",quiet = FALSE, mode = "w")
    fileList <- untar(file.path(destDir,"RSEM_genes_data.tar.gz"),list=TRUE)
    grepSearch = paste("*.*RSEM_genes__data.data.txt",sep="")
    fileList = fileList[grepl(grepSearch,fileList)]
    untar(file.path(destDir,"RSEM_genes_data.tar.gz"),files=fileList)
    fileName<-paste("RSEM_genes_data_",i,".txt",sep="")
    file.rename(from=fileList,to=file.path(destDir,fileName))
    file.remove(file.path(destDir,"RSEM_genes_data.tar.gz"))
    delFodler <- file.path(getwd(),strsplit(fileList,"/")[[1]][1])
    message(delFodler)
    unlink(delFodler, recursive = TRUE)
    
    cat(sprintf("Loading RSEM raw counts\n"))
    
    tmpMat[[i]]<-fread(file.path(destDir,fileName),sep="\t",header=TRUE,stringsAsFactors=FALSE,data.table=FALSE)
    #tmpMat<-read.table(file=paste(destDir,"/","sig_genes.txt",sep=""),header=TRUE,sep="\t")
    
    }
    
    tmpMatCombined<-tmpMat[[1]]
    
    if(length(tmpMat)>1){
     for (j in 2:length(tmpMat)){
        tmpMatCombined<-cbind(tmpMatCombined,tmpMat[[j]][,2:ncol(tmpMat[[j]])])
     }
    }
    
    save(tmpMatCombined,file=file.path(destDir,"RSEM_gene_combined.Rd"))
    
    }else{  
       load(file.path(destDir,"RSEM_gene_combined.Rd"))
    }
    
    return(tmpMatCombined)
    
}
