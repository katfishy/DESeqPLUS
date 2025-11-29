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
            textInput("reference", "Reference", value = "Control"),
            actionButton("run", "Create Plot")
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
      withProgress(message = "Starting DESeqPLUS...", {

        # Step 1: Read files
        counts <- read.csv(input$countFile$datapath, row.names = 1)
        colData <- read.csv(input$coldataFile$datapath, row.names = 1)

        # Make sure they have matching sample names
        # stopifnot(all(colnames(counts) %in% rownames(colData)))

        incProgress(0.3, detail = "Preparing data...")

        # Step 2: Set factors and reference
        design_col <- gsub("~\\s*", "", input$design)
        colData[[design_col]] <- factor(colData[[design_col]])
        colData[[design_col]] <- relevel(colData[[design_col]], ref = input$reference)

        incProgress(0.6, detail = "Creating DESeqDataSet...")

        # Step 3: Create DESeqDataSet
        dds <- DESeqDataSetFromMatrix(countData = counts,
                                      colData = colData,
                                      design = as.formula(input$design))
        incProgress(1, detail = "Done!")

        return(dds)
      })
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
      withProgress(message = "Starting DESeqPLUS...", {
        req(dds_data())
        incProgress(0.5, detail = "Running DESeq2...")
        dds <- DESeq(dds_data())
        res <- results(dds)
        incProgress(0.8, detail = "Creating Plot...")
        DESeqPLUS::volcanoPlot(res)
        incProgress(1, detail = "Done!")
      })
    })
}

# Run the application
shinyApp(ui = ui, server = server)
