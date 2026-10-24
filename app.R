library(shiny)
library(ggplot2)
library(invgamma)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("probability density function plot"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      
        sidebarPanel(        
            selectInput("distribution", label = "distribution",
                        choices = c("normal", "exponential",
                                    "gamma", "inversed gamma")),
            uiOutput("parameters"),
            uiOutput("x_range")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           # textOutput("debug")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$parameters <- renderUI({
    if (input$distribution == "normal") {
      cart = tagList()
      cart[[1]] <- sliderInput("mu",
                  "mu",
                  min = -10,
                  max = 10,
                  value = 0)
      
      cart[[2]] <- sliderInput("sigma",
                  "sigma",
                  min = 0.00001,
                  max = 10,
                  value = 1)
      
      cart
    } else {
    
      if (input$distribution == "exponential") {
        sliderInput("lambda",
                    "rate",
                    min = 0.00001,
                    max = 10,
                    value = 1)
      } else {
      if (input$distribution == "gamma") {
        cart = tagList()
        cart[[1]] <- sliderInput("shape",
                                 "shape",
                                 min = 0.00001,
                                 max = 5,
                                 value = 1)
        
        cart[[2]] <- sliderInput("scale",
                                 "scale",
                                 min = 0.00001,
                                 max = 5,
                                 value = 1)
        
        cart
      } else {
        if (input$distribution == "inversed gamma") {
          cart = tagList()
          cart[[1]] <- sliderInput("shape",
                                   "shape",
                                   min = 0.00001,
                                   max = 5,
                                   value = 1)
          
          cart[[2]] <- sliderInput("scale",
                                   "scale",
                                   min = 0.00001,
                                   max = 5,
                                   value = 1)
          
          cart
      }
      }
    }
  }
  })
  
  output$x_range <- renderUI({
    if (input$distribution == "normal") {
      sliderInput("x_range", "x axis range",
                  value = c(-10, 10), min = -50, max = 50)
    } else {
      if (input$distribution %in% c("exponential", "gamma", "inversed gamma")) {
        sliderInput("x_range", "x axis range",
                    value = c(0, 10), min = 0, max = 50)
      }
      }
  })
    

  
    output$distPlot <- renderPlot({
      if (input$distribution == 'normal'){
        x_grid <- seq(input$x_range[1], input$x_range[2], length.out = 100)
        pdf <- dnorm(x_grid, mean = input$mu, sd = input$sigma)
      }

      if (input$distribution == 'exponential'){
        x_grid <- seq(input$x_range[1], input$x_range[2], length.out = 100)
        pdf <- dexp(x_grid, rate = input$lambda)
      }
      
      if (input$distribution == 'gamma'){
        x_grid <- seq(input$x_range[1], input$x_range[2], length.out = 100)
        pdf <- dgamma(x_grid, shape = input$shape, scale = input$scale)
      }
      
      if (input$distribution == 'inversed gamma'){
        x_grid <- seq(input$x_range[1], input$x_range[2], length.out = 100)
        pdf <- invgamma::dinvgamma(x_grid, shape = input$shape, scale = input$scale)
      }
      
      forplot <- data.frame(x_grid, pdf)
      ggplot(forplot, aes(x_grid, pdf)) +
        geom_line() +
        theme_bw() +
        xlab("x") + ylab("pdf")
    })
    
    # output$debug <- renderText({
    #   paste0(input$x_range3[1])
    # })
}


# Run the application 
shinyApp(ui = ui, server = server)
