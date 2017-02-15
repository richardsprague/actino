# test-phyloseq.r
# test that it's a proper phyloseq object

library(phyloseq)
context("Prove it's a proper phyloseq object")

DATA_DIR <- system.file("extdata", package = "actino") # "../../inst/extdata"

data("kombucha")
p <- kombucha

test_that(paste("Current directory=",getwd()),{
  print("running phyloseq tests now")
  print(getwd())
  #expect_equal(TRUE,FALSE)
})



test_that("Phyloseq object has correct OTU value",{
  expect_equal(as.numeric(otu_table(p)[10,3]),2786) # just an arbitrary element
})

test_that("Phyloseq object has correct number of samples",{
  expect_equal(nsamples(p),19)
})

test_that("create Phyloseq object from a directory of json files",{
  expect_equal(length(just_json_files_in(DATA_DIR)),19)
  pj <- phyloseq_from_JSON_at_rank(just_json_files_in(DATA_DIR),paste0(DATA_DIR,"/kombucha-mapfile.xlsx"))
  expect_equal(nsamples(pj),nsamples(p))
  expect_equivalent(p,pj)

})
