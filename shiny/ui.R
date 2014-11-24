library(shiny)

#load mixed exponential parameter tables
tables <- read.csv("MixExpTables.csv", header = TRUE)

#user interface
shinyUI(fluidPage(
  # Application title
  titlePanel("ILF Sampling Error"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput("table", "Table:", tables$Table.Name),
      numericInput("simulations", "Number of simulations:", 10),
      numericInput("sample.size", "Number of observations per sample:", 50),
      numericInput("limit", "ILF Limit:", 250000),
      submitButton("Update View")),
    mainPanel(
      h3(textOutput("caption")),
      plotOutput("aggPlot"),
      verbatimTextOutput("summary")
    )
  )
))