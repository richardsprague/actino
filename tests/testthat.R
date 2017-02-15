# All test scripts are called from here
# Note that first we read the data that will be used by all test scripts.

library(testthat)
library(actino)

#test_check("actino")

print("running testthat.R")

c <- just_json_files_in("./data/kombucha")

data("kombucha.csv")
data("kombucha")
k = kombucha.csv
k$ssr <- 0000
mapfile <- data.frame(ssr = c(0), Username = c("testPerson"))

p <- kombucha

#print(getwd())
#test_file("./tests/testthat/test-reading.R")
test_dir("./tests/testthat/")

