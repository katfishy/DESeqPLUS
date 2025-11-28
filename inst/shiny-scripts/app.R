library(shiny)
library(DESeq2)
library(ggplot2)

options(shiny.maxRequestSize = 100 * 1024^2)

ui <- fluidPage(
    titlePanel("DESeqPLUS"),

    sidebarLayout(
        sidebarPanel(
            fileInput("countFile", "Upload Count Matrix (CSV)", accept = ".csv"),
            fileInput("coldataFile", "Upload ColData (CSV)", accept = ".csv"),
            textInput("design", "Design Formula", value = "~ condition"),
            textInput("reference", "Reference", value = "reference"),
            actionButton("run", "Run DESeq2")
        ),

        mainPanel(
           tabsetPanel(
             tabPanel("qcBarPlot", plotOutput("qcBarPlot")),
             tabPanel("qcBoxPlot", plotOutput("qcBoxPlot")),
             tabPanel("PCA", plotOutput("pcaPlot")),
             tabPanel("Volcano", plotOutput("volcanoPlot"))
           )
        )
    )
)

server <- function(input, output) {

    dds_data <- eventReactive(input$run, {
      req(input$countFile, input$coldataFile)

      counts <- read.csv(input$countFile$datapath, row.names = 1)
      colData <- read.csv(input$coldataFile$datapath, row.names = 1)

      # Make sure they have matching sample names
      stopifnot(all(colnames(counts) %in% rownames(colData)))

      colData$condition <- factor(colData$condition)
      colData$condition <- relevel(colData$condition, ref = input$reference)

      dds <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = colData,
                                    design = as.formula(input$design))

      return(dds)
    })

    output$qcBarPlot <- renderPlot({
      req(dds_data())
      DESeqPLUS::qcPlot(dds_data())[[1]]
    })

    output$qcBoxPlot <- renderPlot({
      req(dds_data())
      DESeqPLUS::qcPlot(dds_data())[[2]]
    })

    output$pcaPlot <- renderPlot({
      req(dds_data())
      DESeqPLUS::pcaPlot(dds_data())
    })

    output$volcanoPlot <- renderPlot({
      req(dds_data())
      dds <- DESeq(dds_data())
      res <- results(dds)
      DESeqPLUS::volcanoPlot(res)
    })
}

# Run the application
shinyApp(ui = ui, server = server)
