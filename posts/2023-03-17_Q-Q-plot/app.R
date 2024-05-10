library(shiny)
library(ggplot2)
library(stats)

# Define UI
ui <- fluidPage(
  # Use div and style for side-by-side layout
  div(style = "display: flex; justify-content: space-between;",
      div(style = "width: 50%; padding: 10px;", 
          plotOutput("plot1"),
          radioButtons("choice1", "Do you think this sample follows a normal distribution?",
                       choices = c("Yes", "No"), selected = NULL),
          actionButton("submit1", "Submit"),
          textOutput("result1")
      ),
      div(style = "width: 50%; padding: 10px;", 
          plotOutput("plot2"),
          radioButtons("choice2", "Do you think this sample follows a normal distribution?",
                       choices = c("Yes", "No"), selected = NULL),
          actionButton("submit2", "Submit"),
          textOutput("result2")
      )
  )
)

# Define Server Logic
server <- function(input, output, session) {
  # Reactive event for the first plot's submission
  generateData1 <- eventReactive(input$submit1, {
    n <- 30
    dist_number <- floor(runif(1, 0, 3))  # Randomize distribution for fun and independence
    if (dist_number == 0) {
      list(data = rnorm(n), dist = "Normal Distribution")  # Normal distribution
    } else if (dist_number == 1) {
      list(data = rt(n, df = 5), dist = "t Distribution (df=5)")  # t distribution
    } else {
      list(data = runif(n), dist = "Uniform Distribution")  # Uniform distribution
    }
  }, ignoreNULL = FALSE)
  
  # Reactive event for the second plot's submission
  generateData2 <- eventReactive(input$submit2, {
    n <- 30
    dist_number <- floor(runif(1, 0, 3))  # Randomize distribution for fun and independence
    if (dist_number == 0) {
      list(data = rnorm(n), dist = "Normal Distribution")  # Normal distribution
    } else if (dist_number == 1) {
      list(data = rt(n, df = 5), dist = "t Distribution (df=5)")  # t distribution
    } else {
      list(data = runif(n), dist = "Uniform Distribution")  # Uniform distribution
    }
  }, ignoreNULL = FALSE)
  
  # Assign data for plots to avoid unnecessary recalculations
  plotData1 <- reactive({ generateData1()$data })
  plotData2 <- reactive({ generateData2()$data })
  
  # Render Q-Q plot without CI for the first plot
  output$plot1 <- renderPlot({
    data <- data.frame(x = plotData1())
    ggplot(aes(sample = x), data = data) +
      geom_qq_line(col = "steelblue") +
      geom_qq()
  })
  
  # Render Q-Q plot with CI for the second plot
  output$plot2 <- renderPlot({
    data_CI <- data.frame(x = plotData2())
    n <- length(data_CI$x)
    X <- matrix(rnorm(n*1000), ncol = 1000)
    X <- apply(X, 2, sort)
    Q <- apply(X, 1, quantile, probs = c(0.025, 0.975))
    xx <- qnorm(((1:n) - 0.5)/n)
    
    ggplot(aes(sample = x), data = data_CI) +
      geom_qq() +
      geom_ribbon(aes(x = xx, ymin = Q[1, ], ymax = Q[2, ]), alpha = 0.1) +
      geom_line(aes(x = xx, y = xx), col = "steelblue") +
      xlab("Theoretical Quantiles") +
      ylab("Sample Quantiles")
  })
  
  # Results and logic for the first plot
  observeEvent(input$submit1, {
    data_info <- generateData1()
    data <- data.frame(x = data_info$data)
    test_result1 <- shapiro.test(data$x)
    dist_used <- data_info$dist
    response1 <- if (input$choice1 == "Yes") {
      if (test_result1$p.value > 0.05) {
        paste("P-value is high. Data likely follows a normal distribution. True distribution used:", dist_used)
      } else {
        paste("P-value is low. Data likely does not follow a normal distribution. True distribution used:", dist_used)
      }
    } else {
      if (test_result1$p.value <= 0.05) {
        paste("P-value is low. Data likely does not follow a normal distribution. True distribution used:", dist_used)
      } else {
        paste("P-value is high. Data likely follows a normal distribution. True distribution used:", dist_used)
      }
    }
    output$result1 <- renderText(response1)
  })
  
  # Results and logic for the second plot
  observeEvent(input$submit2, {
    data_CI_info <- generateData2()
    data_CI <- data.frame(x = data_CI_info$data)
    test_result2 <- shapiro.test(data_CI$x)
    dist_used_CI <- data_CI_info$dist
    response2 <- if (input$choice2 == "Yes") {
      if (test_result2$p.value > 0.05) {
        paste("P-value is high. Data likely follows a normal distribution. True distribution used:", dist_used_CI)
      } else {
        paste("P-value is low. Data likely does not follow a normal distribution. True distribution used:", dist_used_CI)
      }
    } else {
      if (test_result2$p.value <= 0.05) {
        paste("P-value is low. Data likely does not follow a normal distribution. True distribution used:", dist_used_CI)
      } else {
        paste("P-value is high. Data likely follows a normal distribution. True distribution used:", dist_used_CI)
      }
    }
    output$result2 <- renderText(response2)
  })
}

shinyApp(ui = ui, server = server)

shinyApp(ui = ui, server = server)
