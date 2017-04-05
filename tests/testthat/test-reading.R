# test_reading.R
# tests that prove you can read a file
# assumes you already loaded kombucha.csv and mapfile data

library(actino)
context("Prove you can read data")

data("kombucha.csv")
k <- kombucha.csv

test_that("sample CSV file can be read",{
  expect_equal(exists("kombucha.csv"),TRUE)
  expect_equal(exists("k"),TRUE)
  expect_equal(k[k$tax_name=="Firmicutes",]$count, 47447)
  expect_equal(k[k$tax_name == "root",]$count, 77768)
})

#print(list.files("../../inst/extdata"))

test_that("object can be created from csv file",{
  c <- read.csv(paste0(system.file("extdata", package = "actino"),"/ubiome-export-74607.csv"))
  c$ssr <- 0000
  mapfile <- data.frame(ssr = c(0), Username = c("testPerson"))
  p <- experiment_to_phyloseq(c,mapfile)
  expect_equal(as.numeric(otu_table(p)[11]),1540)
  expect_equal(colnames(tax_table(p)), "Genus")
  expect_equal(rownames(tax_table(p))[1:3],c("g__Acetitomaculum","g__Acidaminococcus","g__Actinomyces"))

})

test_that("can tell which json files are available",{
  expect_equal(length(just_json_files_in(system.file("extdata", package = "actino"))),19)
})

