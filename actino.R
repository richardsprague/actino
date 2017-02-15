# actino.R
# loads the data needed by the app

kombucha.csv <- read.csv("./data/kombucha/ubiome-export-74607.csv")

#devtools::use_data(kombucha.csv, overwrite = TRUE)

#devtools::document()

library(dplyr)
library(tidyr)
library(phyloseq)
library(readxl)
library(actino)
library(jsonlite)

debugSource("R/actinoCSV.r")
debugSource("R/actinoJSON.r")

#kombucha.csv$ssr <- "000"
kombucha.csv$ssr <- 0000
mapfile <- data.frame(ssr = c(0), Username = c("testPerson"))


# k.genus <- kombucha.csv %>% filter(tax_rank=="genus")
# k <- k.genus %>% select(tax_name,ssr,count) %>% spread(tax_name,count)
# k[is.na(k)] <- 0
#
# #k <- kombucha.csv %>% filter(tax_rank=="genus") %>% select(tax_name, ssr, tax_rank, count)  %>% spread(tax_name,count)
#
#
# k <- samples_as_column_dataframe(k.genus)

p <- experiment_to_phyloseq(kombucha.csv,mapfile)

#kombucha.all <- join_all_ubiome_files(just_json_files_in("./data/kombucha"))
kombucha <- phyloseq_from_JSON_at_rank(just_json_files_in("./data/kombucha"),"./data/kombucha/kombucha-mapfile.xlsx")

#devtools::use_data(kombucha, overwrite = TRUE)

auto_test("./R","./tests/testthat")
