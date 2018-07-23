#' @title Select subset of PathwayCommon2 network  
#' 
#' @description 
#' This function provides options to select subset of PathwayCommon2 based on
#' the interaction types or interaction sources.  
#'
#' @details
#' 
#' @param interactionTypes a vector, names of interaction.
#' [ controls-state-change-of,
#'   controls-expression-of,
#'   controls-phosphorylation-of,
#'   controls-transport-of,
#'   catalysis-precedes,
#'   in-complex-with, 
#'   interacts-with  ]
#'
#' @param sources
#'            
#' @return a data frame contains nodes and interactions in four columns
#' [ PARTICIPANT_A, INTERACTION_TYPE, PARTICIPANT_B, INTERACTION_DATA_SOURCE]
#' If INTERACTION_DATA_SOURCE is missing, it will be filled with "NA".  
#'
#' @author Eric Minwei Liu, \email{emliu.research@gmail.com}
#'         
#' @examples 
#' 
#' \dontrun{
#' extendedSifFileNames<-getPC2NetworkName(fileType="EXTENDED_BINARY_SIF")
#' selectedFileName<-extendedSifFileNames[grepl("Reactome",extendedSifFileNames)]
#' extendedSifNetwork <- downloadPc2(selectedFileName=selectedFileName)
#' 
#' names(extendedSifNetwork$edges)
#' 
#' interactionTypes<-unique(extendedSifNetwork$edges$INTERACTION_TYPE)
#' interactionDataSource<-unique(extendedSifNetwork$edges$INTERACTION_DATA_SOURCE)
#' 
#' 
#' 
#' }
#' 
#'
#'            
#' @concept netboxr
#' @export
#' 
preProcessPC2Network <- function(extendedSifNetwork,interactionTypes="EXTENDED_BINARY_SIF",sources=) {
  tmpVector<-extendedSifNetwork$edges$INTERACTION_DATA_SOURCE
  for(i in 1:length(tmpVector)){
    #i<-96
    interactionDataSource[i]<-paste(tmpVector[i][[1]],collapse=";")
  }
  
  tmpDF<-as.data.frame(extendedSifNetwork$edges)
  sifExtra<-tmpDF[,1:3]
  sifExtra$INTERACTION_DATA_SOURCE<-unlist(interactionDataSource)
  
}  