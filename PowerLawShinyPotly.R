#botones
#http://shiny.rstudio.com/gallery/widget-gallery.html

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
#my_palette <- colorRampPalette(c("red","yellow","green","blue"))(n=6)

ui <- fluidPage(
  titlePanel("Uploading Files"),
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose file to upload',
                accept = c(
                  'text/csv',
                  'text/comma-separated-values',
                  'text/tab-separated-values',
                  'text/plain',
                  '.csv',
                  '.tsv'
                )
      ),
      tags$hr(),
      tags$hr(),
      textInput(inputId = "texto", label = "Filtrar por modelo que contenga:"),
      numericInput(inputId = "ancho", label = "Tamaño de banda", value = 3.5), 
      checkboxGroupInput(inputId = "rxn_metab", label = "Escoge rxn/metas o In/Out:", 
                         c("reactions" = "_rxn_", "metabolites" = "_metab_", 
                           "InDegree" = "_In", "OutDegree" = "_Out")),
      radioButtons(inputId = "aic_bic", label = h3("Akaike y/o Bayesian"), 
                   choices = list("Akaike" = "_aic", "Bayesian" = "_bic", "both" = "_aic|_bic"), selected = "_aic|_bic"),
      p('Esta primera version permite controlar',
        'los modelos y el ancho de banda para los graficos',
        a(href = "Javeriana", 'http://www.javeriana.edu.co/')
      )
    ),
    mainPanel(
      tabsetPanel(
        #tabPanel("Tabla", tableOutput('contents')),
        tabPanel("Modelo", tableOutput("texto")),
        tabPanel("Resumen", verbatimTextOutput("summary")),
        tabPanel("graficos", plotlyOutput("graf"))
      )
    )
  )
)

server <- function(input, output) {
  archivo <- reactive({
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    read.csv(inFile$datapath, header = TRUE, sep = " ", quote = '"')
  })
  output$contents <- renderTable({
    if(is.null(archivo())){return(NULL)}
    input$file1
  })
  output$summary <- renderPrint({
    if(is.null(archivo())){return(NULL)}
    tbl_df(archivo()) %>%
      filter(grepl(input$texto, model), grepl(input$rxn_metab, model)) %>%
      summarise(n())
  })
  output$texto <- renderTable({
    tbl_df(archivo()) %>% 
      filter(grepl(input$texto, model))
    
    #print(length(input$rxn_metab))
  })
  dato <- reactive({
    if (length(input$rxn_metab) == 1) {
      input$rxn_metab
    } else if(length(input$rxn_metab) == 2) {
      gsub(" ", "", paste(input$rxn_metab[1], "|", input$rxn_metab[2]))
    } else if (length(input$rxn_metab) == 3) {
      gsub(" ", "", paste(input$rxn_metab[1], "|", input$rxn_metab[2], "|", input$rxn_metab[3]))
    } else if (length(input$rxn_metab) == 4) {
      gsub(" ", "", paste(input$rxn_metab[1], "|", input$rxn_metab[2], "|", input$rxn_metab[3], "|", input$rxn_metab[4]))
    } else {
      "_"
    }
    
  })
  output$graf <- renderPlotly({
    #print(dato())
    grafo <- archivo() %>%
      dplyr::select(ends_with(".y"), model) %>%
      filter(grepl(input$texto, model), grepl(dato(), model), grepl(input$aic_bic, model)) %>% 
      gather(distribucion, valor, PowerLaw.y:LogLogistic.y) %>%
      mutate(in_out_total = ifelse(grepl('_Total', model) == TRUE, 'Total', ifelse(grepl('_In', model) == TRUE, 'In', 'Out'))) %>%
      ggplot(aes(x = model, y = distribucion, fill = valor)) + geom_tile(color = "black", width = input$ancho) + facet_grid(~in_out_total)
    ggplotly(grafo)
    #ggplot(aes(x = model, y = distribucion, fill = factor(valor))) + geom_tile( width = input$ancho) + facet_grid(~in_out_total) + scale_fill_manual(values = my_palette)
  })
}

shinyApp(ui = ui, server = server)