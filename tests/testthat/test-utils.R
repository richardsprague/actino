# test-utils.r
# test utilities for manipulating uBiome files

context("Test utilities for manipulating uBiome files")

#library(actino)
DATA_DIR <- system.file("extdata", package = "actino") # "../../inst/extdata"


a <- join_all_ubiome_files_full(just_json_files_in(DATA_DIR)[1:2])
a[is.na(a)] <- 0

ssrs<-sapply(strsplit(names(a)[c(-1)],"\\$"),function(x) as.numeric(x[2]))
p <- kombucha

test_that("It's a dataframe with correct SSRs",{
  expect_equal(is.data.frame(a),TRUE) # just an arbitrary element
  expect_equal(ssrs,c(73813,73837))
  #expect_equal(a$tax_name,rep("Rkoot",237))
  #print(a[3,1])
  expect_equal(a[3,1],"Root;n__cellular organisms;k__Bacteria;l__Bacteroidetes/Chlorobi group;p__Bacteroidetes;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae;g__Flavobacterium")
  expect_equal(a[15,1],"Root;n__cellular organisms;k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia")

})

test_that("matrix at rank shows correct items",{
  m <- matrix_at_rank(a, rank = "genus")
  expect_equal(nrow(m),83)
  expect_equal(m[6,2],4300)
  expect_equal(matrix_at_rank(a, rank = "Root")[1,1],67550)
})
