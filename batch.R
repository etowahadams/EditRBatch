library(gamlss)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)

# Parameters
key.filename <- "230105_ExampleData_Key.xlsx"
seqfile.dir <- "ExampleDataFiles"

data.key <- read_excel(key.filename)
data.key %<>% fill(gRNA, "guide sequence", "Reverse Y/N")
# data.key <- data.key[10,]

data.key$results <- apply(data.key, 1, function (x) {
  tryCatch({
    CalcEditRBatch(x[[1]], seqfile.dir, x[[3]], x[[4]])
  }, error=function(cond) {
    return(paste("Error:", cond))
  })
})