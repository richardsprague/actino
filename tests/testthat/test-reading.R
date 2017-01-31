# testReading.R
# tests that prove you can read a file

library(actino)
context("Prove you can read data")

data("kombucha.csv")
k = kombucha.csv
k$ssr <- 0000

test_that("sample CSV file can be read",{
  expect_equal(exists("kombucha.csv"),TRUE)
  expect_equal(k[k$tax_name=="Firmicutes",]$count, 47447)
})



#k <- read.csv("./data/kombucha/ubiome-export-74607.csv")

