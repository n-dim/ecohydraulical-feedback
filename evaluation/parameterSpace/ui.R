#---- load parameters ----
load(file="../../example simulation run/exampleParameters.RData")

parlist$pa <- c(100,200,300)
parlist$roughness <- 10^(-2:2)

values=NA

#---- calculate parameter space----

Dimnames <- parlist
Dims <- mapply(length, Dimnames)
totalDims <- prod(Dims)
parameterSpace <- array(1:prod(Dims), dim=Dims, Dimnames)

library(R.utils)
selectedSims <- extract(parameterSpace, run="T",  drop=T)


# Define UI for dataset viewer application
shinyUI(pageWithSidebar(
  
  # Application title
  headerPanel("DWD Daten"),
  
  # Sidebar with controls to provide a caption, select a dataset, and 
  # specify the number of observations to view. Note that changes made
  # to the caption in the textInput control are updated in the output
  # area immediately as you type
  sidebarPanel(
    selectInput("sel1", "select parameter:", 
                choices = names(Dimnames)),
    uiOutput("selection2"),

  
    numericInput("time", "Timestep:", 10)
  ),
  
  
  # Show the caption, a summary of the dataset and an HTML table with
  # the requested number of observations
  mainPanel(
    h3(textOutput("caption")), 
    
    #verbatimTextOutput("summary"), 
    
    #tableOutput("view")
    
    mainPanel(
      verbatimTextOutput("parameterSpace")
    )
  )
))
