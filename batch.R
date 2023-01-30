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

CalcEditRBatch <- function(file.id, file.dir, guideseq, reverseYN) {
  p.val.cutoff <- 0.01
  default.trim <- FALSE
  is.reverse <- (reverseYN == "Y")
  
  # get filename
  dir.files <- list.files(file.dir)
  file.matches <- grep(paste0("^", file.id, "_"), dir.files)
  if (length(file.matches) != 1) {
    stop(paste("The provided file ID", file.id, "matched multiple files in the directory or could not be found"))
  }
  filename <- paste0(file.dir, "/", dir.files[[file.matches[[1]]]])
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
  if (is.reverse) {
    results$table %<>% filter(guide.seq == "G") %>% filter(focal.base == "A")
  } else {
    results$table %<>% filter(guide.seq == "C") %>% filter(focal.base == "T")
  }
  results
}

data.key$results <- apply(data.key, 1, function (x) {
  tryCatch({
    CalcEditRBatch(x[[1]], seqfile.dir, x[[3]], x[[4]])
  }, error=function(cond) {
    return(paste("Error:", cond))
  })
})


out <- apply(data.key, 1, function(x) {
  result_table <- x[[6]][[1]]
  if (length(result_table) > 2) {
    result_table %<>% arrange(guide.position)
    result_table[["5 trim"]] <- x[[6]][[3]][[1]]
    result_table[["3 trim"]] <- x[[6]][[3]][[2]]
    result_table[["Sample ID"]] <- x[[1]]
    result_table[["gRNA"]] <- x[[2]]
    result_table[["guide sequence"]] <- x[[3]]
  } else {
    result_table[[12]] <- x[[1]]
    result_table[[13]] <- x[[2]]
    result_table[[14]] <- x[[3]]
  }
  result_table
})

out <- do.call(rbind, out)
# redo the order so it makes sense 
out <- out[, c(14, 13, 12, 5, 8, 1, 2, 3, 4, 6, 7, 9, 10, 11)] 



