#' Making oncoprint
#' 
#' @param alterationMat data frame of alteration
#' @param cancerType TCGA cancer study name
#' @param moduleNum module number
#' @param workDir Path to the location of file save
#' @param width Default is 15 inch
#' @param height Default is 10 inch
#' @param verbose Default is TRUE
#' 
#' @return null
#'   
#' @examples
#' date<-"2015_08_21"
#' cancerType<-"KIRC"
#' worDir<-getwd()
#' #mutationMAF<-getMutationMAF(date,cancerType,workDir)
#' 
#' #alterationMatSingle<-unionMultipleFrames(alterationMatMultiple)
#' 
#' @concept netboxr
#' @export 
#' @importFrom ComplexHeatmap oncoPrint
#' @import GetoptLong
makeOncoprint<-function(alterationMat,cancerType,moduleNum,workDir=NULL,width=15,height=10,verbose=TRUE){

#filePath<-"~/work/Ekta_lab/COADREAD_project"
#fileName<-"tawny2.txt"
#alterationMat<-read.table(paste(filePath,fileName,sep="/"),sep="\t",header=TRUE,
#                          stringsAsFactors=FALSE,fill=TRUE,quote="\"", comment='',na.strings=c("","NA"))

mat2<-alterationMat
#mat<-mat[which(rownames(mat) %in% selectedGene),]

mat2[is.na(mat2)]<-""
mat2<-as.matrix(mat2)

mat<-alterationMat

#####

matOriginal <- alterationMat

# remove sample column does not have any alterations
#mat = mat[, !apply(mat, 2, function(x) all(grepl("^\\s*$", x)))]
matReduced = matOriginal[, !apply(matOriginal, 2, function(x) all(is.na(x)))]

# alteraed ratio is ( size of sample with alterations ) / ( total sample size )
altered <- ncol(matReduced)/ncol(matOriginal)


#####


alterationScore<-{}
alterType<-c("MUT","AMP","HOMDEL","EpiSilenced")
alterationScore[[alterType[1]]]<-1
alterationScore[[alterType[2]]]<-3
alterationScore[[alterType[3]]]<-2
alterationScore[[alterType[4]]]<-4

vectorTmp<-unlist(mat)
vectorTmp<-as.character(vectorTmp)
for(i in 1:length(vectorTmp)){
  
  content<-strsplit(as.character(vectorTmp[i]),";")[[1]]
  
  score<-0
  for(j in 1:length(content)){
    if( content[j] %in% alterType){
      score<-score + alterationScore[[content[j]]]
    }
  }
  vectorTmp[i]<-score
  
}

matTmp<-matrix(as.numeric(vectorTmp),nrow=nrow(mat))

memoSort<-function(matTmp){
  
  #matTmp<-matTmp>0
  matTmp2<-(matTmp>0)
  
  # sorting on the sample size of alterations per gene
  #geneOrder <- sort(rowSums(mat), decreasing=TRUE, index.return=TRUE)$ix;
  geneOrder <- sort(rowSums(matTmp2), decreasing=TRUE, index.return=TRUE)$ix
  
  
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]>0) {
        #score <- score + 2^(length(x)-i)
        #score <- score + 5^(length(x)-i)
        #score <- score + 5^(length(x)-i) + x[i]
        #score <- score + (x[i]+1)*10^(length(x)-i)
        score <- score + (x[i]+1)*5^(length(x)-i) 
        #score <- score + ( (x[i]+1) /length(x) ) * ( 2^(length(x)-i) ) 
        #score <- score + (length(x)-i)^100 + (length(x)-i)^(x[i]+1)
      }
    }
    return(score);
  }
  
  scores <- apply(matTmp[geneOrder, ], 2, scoreCol)
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix
  memoSortOrder<-list(geneOrder=geneOrder,sampleOrder=sampleOrder,scores=scores)
  
  return(memoSortOrder)
}

memoSortOrder<-memoSort(matTmp)

#####

tmpVector<-unlist(mat)
alterTypes<-names(table(as.character(tmpVector)))
alterTypes<-unique(unlist(strsplit(alterTypes,";")))

alterFactor<-c(1,2,3,4)
names(alterFactor)<-c("MUT","AMP","HOMDEL","EPISILENCED")

alterNumber<-sum(alterFactor[alterTypes])

backgroundFun<-function(x, y, w, h) {
  grid.rect(x, y, w-unit(0.0, "mm"), h-unit(0.0, "mm"), gp = gpar(fill = "white", col = "#CCCCCC"))
}

epiSilencedFun<-function(x, y, w, h) {
  grid.rect(x, y, w-unit(0.0, "mm"), h-unit(0.0, "mm"), gp = gpar(fill = "black", col = "#CCCCCC"))
} 

homdelFun<-function(x, y, w, h) {
  grid.rect(x, y, w-unit(0.0, "mm"), h-unit(0.0, "mm"), gp = gpar(fill = "blue", col = "#CCCCCC"))
}  

ampFun = function(x, y, w, h) {
  grid.rect(x, y, w-unit(0.0, "mm"), h-unit(0.0, "mm"), gp = gpar(fill = "red", col = "#CCCCCC"))
}

mutFun = function(x, y, w, h) {
  grid.rect(x, y, w-unit(0.0, "mm"), h*0.33, gp = gpar(fill = "#008000", col = "#CCCCCC"))
}


if(length(alterTypes)==1){
  alterFunctionList<-switch(alterTypes,
                            "MUT"=list(background = backgroundFun, MUT=mutFun),
                            "AMP"=list(background = backgroundFun, AMP=ampFun),
                            "HOMDEL"=list(background = backgroundFun, HOMDEL=homdelFun),
                            "EPISILENCED"= list(background = backgroundFun, EPISILENCED=epiSilencedFun)   
  )
}

if(length(alterTypes)==2){
  
  alterName<-{}
  alterFunctionList<-{}
  
  if(alterNumber==3){alterName<-"MA"}
  if(alterNumber==4){alterName<-"MH"}
  if(alterNumber==5 & ("MUT" %in% alterTypes)){alterName<-"ME"}
  if(alterNumber==5 & ("AMP" %in% alterTypes)){alterName<-"AH"}
  if(alterNumber==6){alterName<-"AE"}
  if(alterNumber==7){alterName<-"HE"}
  
  
  alterFunctionList<-switch(alterName,
                            "MA"=list(background = backgroundFun, 
                                      MUT=mutFun,
                                      AMP=ampFun),
                            "ME"=list(background = backgroundFun,
                                      MUT=mutFun,
                                      EPISILENCED=epiSilencedFun),
                            "MH"=list(background = backgroundFun, 
                                      MUT=mutFun,
                                      HOMDEL=homdelFun),
                            "AH"= list(background = backgroundFun, 
                                       AMP=ampFun,
                                       HOMDEL=homdelFun),
                            "AE"=list(background = backgroundFun, 
                                      AMP=ampFun,
                                      EPISILENCED=epiSilencedFun),
                            "HE"= list(background = backgroundFun, 
                                       HOMDEL=homdelFun,
                                       EPISILENCED=epiSilencedFun)
  )
}

if(length(alterTypes)==3){
  alterName<-{}
  alterFunctionList<-{}
  
  if(alterNumber==6){alterName<-"MAH"}
  if(alterNumber==7){alterName<-"MAE"}
  if(alterNumber==8){alterName<-"MHE"}
  if(alterNumber==9){alterName<-"AHE"}
  
  alterFunctionList<-switch(alterName,
                            "MAH"=list(background = backgroundFun, 
                                       MUT=mutFun,
                                       AMP=ampFun,
                                       HOMDEL=homdelFun),
                            "MAE"=list(background = backgroundFun,
                                       MUT=mutFun,
                                       AMP=ampFun,
                                       EPISILENCED=epiSilencedFun),
                            "MHE"=list(background = backgroundFun, 
                                       MUT=mutFun,
                                       HOMDEL=homdelFun,
                                       EPISILENCED=epiSilencedFun),
                            "AHE"= list(background = backgroundFun, 
                                        AMP=ampFun,
                                        HOMDEL=homdelFun,
                                        EPISILENCED=epiSilencedFun)   
  )
}

if(length(alterTypes)==4){
  alterFunctionList<-list(background = backgroundFun, 
                          MUT=mutFun,
                          AMP=ampFun,
                          HOMDEL=homdelFun
  )
  
}

#####

alter_fun_list = alterFunctionList

col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue","EpiSilenced"="black")

dataDir<-workDir
filePath<-file.path(dataDir,cancerType,"oncoprint")

if( !file.exists(filePath) ){
  dir.create(filePath,recursive=TRUE)
}


fileName<-paste(cancerType,"_oncoprint_module_",moduleNum,".pdf",sep="")
fileName<-file.path(filePath,fileName)
pdf(fileName, width = width, height = height)
ht_list<-oncoPrint(mat2, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun_list = alter_fun_list, col = col,
          row_order=memoSortOrder$geneOrder,column_order=memoSortOrder$sampleOrder,
          remove_empty_columns = TRUE,
          column_title = qq("OncoPrint for TCGA @{cancerType} module @{moduleNum}\nAltered in @{ncol(matReduced)}/@{ncol(matOriginal)} (@{round(altered*100)}% of cases)"),
          #column_title_gp = gpar(fontsize = 12, fontface = "bold")),
          heatmap_legend_param = list(title = "Alterations", at = c("MUT", "AMP", "HOMDEL","EpiSilenced"), 
                                      labels = c("MUT", "AMP", "HOMDEL","EpiSilenced")))
draw(ht_list)

dev.off()

}


