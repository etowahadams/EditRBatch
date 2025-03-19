ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Choose ab1 files", multiple=TRUE),
      fileInput("data.key", "Choose metadata file (format: 'Sample ID,	gRNA	guide sequence,	Reverse Y/N
')", multiple=TRUE)
    ),
    mainPanel(
      tableOutput("contents")
    )
  )
)
