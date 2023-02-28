library(readxl)

server <- function(input, output) {
  output$contents <- renderTable({
    filenames <- input$files$name
    datapaths <- input$files$datapath
    keyfile <- input$data.key$datapath
    
    if (!is.null(filenames) & !is.null(keyfile)) {
      
      data.key <- read_excel(keyfile)
      data.key %<>% fill(gRNA, "guide sequence", "Reverse Y/N")
      
      # for each row in the data key, apply the CalcEditRBatch function 
      n <- nrow(data.key)
      all_data <- data.frame(matrix(ncol = 16, nrow = 0))
      output_cols <- c("sample.id", "gRNA.id", "gRNA.seq","reverseYN", "guide.position","perc","pval","Tot.area","base.call","index","guide.seq","focal.base","area", "trim5", "trim3", "error")
      colnames(all_data) <- output_cols
      all_data$gRNA.seq <- as.character(all_data$gRNA.seq)
      all_data$error <- as.character(all_data$error)
      
      withProgress(message = 'Calculating values', value = 0, {
        for (i in 1:n) {
          incProgress(1/n, detail = paste("Analyzing file", i))
          
          sample_info <- data.key[i,]  # row from data key  
          
          sample.id <- sample_info[[1]]
          gRNA.id <- sample_info[[2]]
          gRNA.seq <- sample_info[[3]]
          reverseYN <- sample_info[[4]]
          
          sample_results <- tryCatch({
            result <- CalcEditRBatch(sample.id, filenames, datapaths, gRNA.seq, reverseYN)
            result_table <- result[[1]]
            trim <- result[[3]]
            
            result_table %<>% arrange(guide.position)
            result_table[["trim5"]] <- trim[[1]]
            result_table[["trim3"]] <- trim[[2]]
            # add in the data from the data key
            result_table[["sample.id"]] <- sample.id
            result_table[["gRNA.id"]] <- gRNA.id
            result_table[["gRNA.seq"]] <- gRNA.seq
            result_table[["reverseYN"]] <- reverseYN
            result_table[["error"]] <- ""
            # append to main 
            # all_data <- rbind(all_data, result_table)
            result_table
          }, error=function(cond) {
            new_row <- list(sample.id=sample.id, gRNA.id=gRNA.id, gRNA.seq=gRNA.seq,reverseYN=reverseYN,guide.position=NA,perc=NA,pval=NA,Tot.area=NA,base.call=NA,index=NA,guide.seq=NA,focal.base=NA,area=NA, trim5=NA, trim3=NA, error=cond$message)
            return(new_row)
          })
          # append the sample results to the result
          all_data <- rbind(all_data, sample_results)
        }
      })
      
      # reorder columns
      all_data <- all_data[,output_cols]
      return(all_data)
    } else {
      return("Please select files")
    }
    
  })
}

# Run the application 
# shinyApp(ui = ui, server = server)
