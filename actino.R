# actino.R
# loads the data needed by the app
# Also use this as the starting point if you want to change something
# source this file (beware of what you're uncommenting) and then run the autotest at the end.

library(dplyr)
library(tidyr)
library(phyloseq)
#library(readxl)
#library(actino)
#library(jsonlite)
library(testthat)

# if you want to hand-assemble the functions from the package:
# debugSource("R/actinoCSV.r")
# debugSource("R/actinoJSON.r")

DATA_DIR <- system.file("extdata", package = "actino") # "../../inst/extdata"

kombucha.csv <- read.csv(paste0(DATA_DIR,"/ubiome-export-74607.csv"))

# Uncomment the following lines to save the data to the package
#devtools::use_data(kombucha.csv, overwrite = TRUE)

# Uncomment this line to generate all the documentation
devtools::document()


# Next, assemble a dataframe and associated mapfile and then create a new phyloseq object
kombucha.csv$ssr <- 0000
mapfile <- data.frame(ssr = c(0), Username = c("testPerson"))

p <- experiment_to_phyloseq(kombucha.csv,mapfile)

kombucha.genus <- phyloseq_from_JSON_at_rank(just_json_files_in(DATA_DIR),paste0(DATA_DIR,"/kombucha-mapfile.xlsx"))

kombucha.all_ranks <- phyloseq_from_JSON(just_json_files_in(DATA_DIR),paste0(DATA_DIR,"/kombucha-mapfile.xlsx"))


devtools::use_data(kombucha.genus,kombucha.all_ranks,kombucha.csv, overwrite = TRUE)

# this line is useful for interactive testing:

auto_test("./R","./tests/testthat")
