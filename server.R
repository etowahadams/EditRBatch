library(readxl)

server <- function(input, output) {
  output$contents <- renderTable({
    filenames <- input$files$name
    datapaths <- input$files$datapath
    keyfile <- input$data.key$datapath
    
    if (!is.null(filenames) & !is.null(keyfile)) {
      
      data.key <- read_excel(keyfile)
      data.key %<>% fill(gRNA, "guide sequence", "Reverse Y/N")
      
      n <- nrow(data.key)
      withProgress(message = 'Calculating values', value = 0, {
        data.key$results <- sapply(1:n, function (i) {
          x <- data.key[i,]
          incProgress(1/n, detail = paste("Doing part", i))
          
          tryCatch({
            CalcEditRBatch(x[[1]], filenames, datapaths, x[[3]], x[[4]])
          }, error=function(cond) {
            return(paste("Error:", cond))
          })
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
          result_table[[12]] <- x[[1]] # sample ID
          result_table[[13]] <- x[[2]] # gRNA
          result_table[[14]] <- x[[3]] # guide sequence 
        }
        result_table
      })
      
      out <- do.call(rbind, out)
      # redo the order so it makes sense 
      out <- out[, c(14, 13, 12, 5, 8, 1, 2, 3, 4, 6, 7, 9, 10, 11)] 
      
    } else {
      return("Please select files")
    }
    
  })
}

# Run the application 
# shinyApp(ui = ui, server = server)
