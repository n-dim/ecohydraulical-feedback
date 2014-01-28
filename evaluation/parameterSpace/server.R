library(shiny)
library(datasets)
library(R.utils)

values=NULL

# Define server logic required to summarize and view the selected dataset
shinyServer(function(input, output) {
  
  
  output$caption <- renderText({
    "test"
  })
  output$selection2 <- renderUI({
    values <- c(Dimnames[[input$sel1]], "x-axis", "y-axis")
    selectInput("sel2", "select value:", 
                choices = values)
  })
  output$parameterSpace <- renderPrint({
    extract(parameterSpace, indices=input$sel2, dims=input$sel1,  drop=T)
  })
  

  




})
