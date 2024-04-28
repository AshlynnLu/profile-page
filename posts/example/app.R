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
  data <- rnorm(100)  # You can change this sample data
  
  # Render Q-Q plot
  output$qqPlot <- renderPlot({
    qqnorm(data)
    qqline(data, col = "steelblue")
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
