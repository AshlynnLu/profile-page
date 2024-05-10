library(shiny)
library(ggplot2)
library(stats)

ui <- fluidPage(
  titlePanel("Interactive Q-Q Plot Analysis"),
  sidebarLayout(
    sidebarPanel(
      radioButtons("choice", "Do you think this sample follows a normal distribution?",
                   choices = c("Yes", "No"), selected = NULL),
      actionButton("submit", "Submit")
    ),
    mainPanel(
      plotOutput("qqPlot"),
      textOutput("result")
    )
  )
)

server <- function(input, output, session) {
  # Sample data
  eg_norm_CI <- data.frame(x = rnorm(n)) # You can change this sample data
  
  X <- matrix(rnorm(n*1000), nr = n)
  X <- apply(X, 2, sort)
  Q <- apply(X, 1, quantile, prob = c(0.025, 0.975))
  xx <- qnorm(((1:n)-0.5)/n)
  
  # Render Q-Q plot
  output$qqPlot <- renderPlot({
    ggplot(aes(sample = x), data = eg_norm_CI) +
    geom_qq() +
    geom_ribbon(aes(x = xx, ymin = Q[1, ], ymax = Q[2, ]), alpha = 0.1) +
    geom_line(aes(x = xx, y = xx)) +
    xlab("theoretical") +
    ylab("sample")
  })
  
  observeEvent(input$submit, {
    # Perform Shapiro-Wilk test when user submits their choice
    test_result <- shapiro.test(data)
    
    # Check user's answer and display appropriate feedback
    if (!is.null(input$choice)) {
      response <- if (input$choice == "Yes") {
        if (test_result$p.value > 0.05) {
          "Correct, p-value is high. Data likely follows a normal distribution."
        } else {
          "Incorrect, p-value is low. Data likely does not follow a normal distribution."
        }
      } else {
        if (test_result$p.value <= 0.05) {
          "Correct, p-value is low. Data likely does not follow a normal distribution."
        } else {
          "Incorrect, p-value is high. Data likely follows a normal distribution."
        }
      }
      output$result <- renderText(response)
    }
  })
}

shinyApp(ui = ui, server = server)
