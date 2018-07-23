context("Testing getCancerTypeSummary()")


test_that("TCGA data summary tests", {
  
  ## Simple test
  date<-"2014_10_17"
  cancerTypeSummary<-getCancerTypeSummary(date=date,dataType="all")
  cancerType<-"BRCA"
  mRNASeq<-1091
  query<-cancerTypeSummary[cancerTypeSummary$Cohort %in% cancerType,]$mRNASeq

  expect_equal(query,mRNASeq)
  
  
  
}

)