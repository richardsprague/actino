# test-phyloseq.r
# test that it's a proper phyloseq object

#library(actino)
context("Prove it's a proper phyloseq object")

data("kombucha")
p <- kombucha
test_that("Phyloseq object has correct OTU value",{
  expect_equal(as.numeric(otu_table(p)[10,3]),2786) # just an arbitrary element
})

test_that("Phyloseq object has correct number of samples",{
  expect_equal(nsamples(p),19)
})

test_that("create Phyloseq object from a directory of json files",{
  data_dir <- "../../data/kombucha"
  expect_equal(length(just_json_files_in(data_dir)),19)
  phyloseq_from_JSON_at_rank(just_json_files_in(data_dir),paste0(data_dir,"/kombucha-mapfile.xlsx"))

})
