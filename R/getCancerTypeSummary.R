#' Retrieve TCGA tumor type data summary
#' 
#' @param date Default is "latest" 
#' @param dataType Default is "TP" for primary tumor. [ TP,NT,all,complete ]
#' 
#' @return a data frame of data summary
#'   
#' @examples
#' date<-getRunDates(latest=TRUE)
#' cancerTypeSummary<-getCancerTypeSummary(date=date,dataType="all")
#' cancerTypeSummary
#' 
#' @concept netboxr
#' @export
#' @import XML
getCancerTypeSummary<-function(date="latest",dataType="TP"){

    #library("XML")
    #date<-"2014_12_06"
    if( date %in% "latest"){
      date<-getRunDates(latest=TRUE)
    }
    
    cat(sprintf("Retreive data from %s\n",date))
    
    url<-"http://gdac.broadinstitute.org/runs"
    url<-paste(url,"/stddata__",date,sep="")
    url<-paste(url,"/samples_report",sep="")
    doc<-htmlTreeParse(url,useInternalNodes=TRUE)

    keyWord = paste("","sample_counts",sep="")
    keyWord = paste("//a[contains(@href, '",keyWord,"')]",sep="")
    plinks1 = xpathSApply(doc, keyWord, xmlAttrs)
    data_fn<-plinks1[2]
    
    url<-"http://gdac.broadinstitute.org/runs"
    url<-paste(url,"/stddata__",date,sep="")
    url<-paste(url,"/samples_report",sep="")
    url<-paste(url,"/",data_fn,sep="")
    
    data<-read.table(url,
                     header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

    data<-data[ !data$Cohort %in% "Totals",]
    
    if(dataType %in% "all"){      
       detailedSummaryLine<-grep("-",data$Cohort)
       keep<-setdiff(1:nrow(data),detailedSummaryLine)
       return(data[keep,]) }
    
    if(dataType %in% "NT"){
       keep<-grep(dataType,data$Cohort)
       return(data[keep,])
    }
    
    if(dataType %in% "TP"){
      keep<-grep(dataType,data$Cohort)
      return(data[keep,])
    }
    
    if(dataType %in% "complete"){
      return(data)
    }
    
}





