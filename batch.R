library(gamlss)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)

# Parameters
key.filename <- "230105_ExampleData_Key.xlsx"
seqfile.dir <- "ExampleDataFiles"



data.key <- read_excel(key.filename)
data.key %<>% fill(gRNA, "guide sequence", "Reverse Y/N")

CalcEditRBatch <- function(file.id, file.dir, guideseq, reverseYN) {
  p.val.cutoff = 0.01
  default.trim = FALSE
  is.reverse = if (reverseYN == "Y") TRUE else FALSE
  
  # get filename
  dir.files <- list.files(file.dir)
  file.matches <- grep(paste0("^", file.id, "_"), dir.files)
  if (length(file.matches) != 1) {
    stop(paste("The provided file ID", file.id, "matched multiple files in the directory or could not be found"))
  }
  filename <- paste0(file.dir, "/", dir.files[[file.matches[[1]]]])
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
}

data.key$results <- apply(data.key, 1, function (x) {
  tryCatch({
    CalcEditRBatch(x[[1]], seqfile.dir, x[[3]], x[[4]])
  }, error=function(cond) {
    return(paste("Error:", cond))
  })
})