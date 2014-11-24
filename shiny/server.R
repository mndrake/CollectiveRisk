library(shiny)

#load mixed exponential parameter tables
tables <- read.csv("MixExpTables.csv", header = TRUE)
weights <-read.csv("MixExpWeights.csv", header = FALSE, stringsAsFactors = TRUE)
means <- read.csv("MixExpMeans.csv", header = FALSE, stringsAsFactors = TRUE)

#load libraries
source("CollectiveRiskS3.R")

#ILF sample plot function
sample.data <- function(idx, simulations, sample.size, limit) {
  d <- MixedExponential(weights = as.vector(weights[1:12,idx]),
                        means = as.vector(means[1:12,idx]))
  set.seed(1)
  ILF.samples <- function(simulations, sample.size, limit){
    p.max <- cdf(d,limit)
    sample <- function(){
      samples <- sapply(runif(sample.size), 
                        function(p) {
                          if (p>=p.max){limit}
                          else {icdf(d,p)}})
      mean(samples) / mean(sapply(samples, function(x){min(x,100000)}))
    }  
    replicate(simulations, sample(), simplify = "array")
  }
  
  ILF <- ILF.samples(simulations, sample.size, limit)
  
  ILF
}

shinyServer(function(input, output) {
        
 datasetInput <- function() {
   idx <- which(tables$Table.Name==input$table)
   sample.data(idx, input$simulations, input$sample.size, input$limit)
 }
  
    output$caption <- renderText(input$table)
    
    output$summary <- renderPrint({
        dataset <- datasetInput()
        print(summary(dataset))
        print('quantiles')
        print(quantile(dataset, c(.05,.95)))
        print('% of mean for quantiles')
        m <- mean(dataset)
        p05 <- quantile(dataset, 0.05)
        p95 <- quantile(dataset, 0.95)
        print(c(p05/m-1, p95/m-1))
        print('CoV (want < 0.0607)')
        sd(dataset) / m
      })   
    
   output$aggPlot <- renderPlot({
     d <- density(datasetInput())
     plot(d, main = "Sampling Error")
     polygon(d, col="red", border="blue")
   })
})