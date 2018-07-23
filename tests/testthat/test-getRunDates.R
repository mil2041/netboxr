context("Testing getRunDates()")


test_that("TCGA analysis run date tests", {
  
  ## Simple test
  date <- getRunDates(latest=FALSE)
  
  expect_equal(sum(date %in% "2014_10_17"),1)
  
  
  
}

)