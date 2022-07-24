
install.packages(shiny)
library(shiny)

# Define UI for slider demo app ----
ui <- fluidPage(

  # App title ----
  titlePanel("What do fold differences in absolute abundance look like at various scales?"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar to demonstrate various slider options ----
    sidebarPanel(

      # Input: Simple integer interval ----
      sliderInput("integer_one_new", "Log fold difference (cases/controls) in taxon 1 absolute abundance:",
                  min = -1, max = 1,
                  value = 0,
                  step = .1),
    sliderInput("integer2", "Log fold difference (cases/controls) in taxon 2 absolute abundance:",
                min = -1, max = 1,
                value = 0,
                step = .1),
  sliderInput("integer_three_new", "Log fold difference (cases/controls) in taxon 3 absolute abundance:",
              min = -1, max = 1,
              value = 0,
              step = .1)),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Table summarizing the values entered ----
      # plotOutput(outputId = "distPlot")

      splitLayout(cellWidths = c("33%", "33%","33%"), plotOutput("fold_differences"),plotOutput("distPlot0"), plotOutput("distPlot"))

      )

    )
  )
# )

# Define server logic for slider examples ----
server <- function(input, output) {

  # Reactive expression to create data frame of all input values ----
  sliderValues <- reactive({

    data.frame(
      Name = c("integer_three_new"),
      Value = as.character(c(input$integer)),
      stringsAsFactors = FALSE)

  })

  # Show the values in an HTML table ----
  output$fold_differences <-
    renderPlot({
      plot(x = 1:3,
        y = c(input$integer_one_new,input$integer2,input$integer_three_new),
        xlab = "Taxon",cex = 1.5,
        cex.lab = 1.5,
        ylab = "",
        col = c("#66C2A5","#FC8D62","#8DA0CB"),
        pch = 19,
        ylim = c(-1,1))
      title(ylab="Log Fold Difference (Cases/Controls)", line=2,
            cex.lab=1.5)
    })
  output$distPlot0 <- renderPlot({

    g1 <- rep(6.3*1e6,3)
    g2 <- g1
    g2[3] <- g2[3]*exp(input$integer_three_new)
    g2[2] <- g2[2]*exp(input$integer2)
    g2[1] <- g2[1]*exp(input$integer_one_new)

    # g2 <- g2/sum(g2)

    together <- as.matrix(cbind(g1,g2))
    rownames(together) <- c("Taxon 1","Taxon 2", "Taxon 3")
    colnames(together) <- c("Controls","Cases")

    barplot(together, border="white", xlab="Group",
            col = c("#66C2A5","#FC8D62","#8DA0CB"),
            # ylab = expression(paste("Cell Concentrations (copies/",mu,"L)",
            # collapse = "")),
            ylim = c(0, 2e7),
            cex.lab = 1.5,
    beside = TRUE)
title(ylab=expression(paste("Cell Concentrations (copies/",mu,"L)")), line=2,
                            cex.lab=1.5)

    text(1.5,2.5e6,
         "Taxon 1",cex = 1.1,srt = 90)
    text(2.5,2.5e6,
         "Taxon 2",cex = 1.1,srt = 90)

    text(3.5,2.5e6,
         "Taxon 3",cex = 1.1,srt = 90)

  })
  output$distPlot <- renderPlot({

    g1 <- rep(1,3)/3
    g2 <- g1
    g2[3] <- g2[3]*exp(input$integer_three_new)
    g2[2] <- g2[2]*exp(input$integer2)
    g2[1] <- g2[1]*exp(input$integer_one_new)

    g2 <- g2/sum(g2)

    together <- as.matrix(cbind(g1,g2))
    rownames(together) <- c("Taxon 1","Taxon 2", "Taxon 3")
    colnames(together) <- c("Controls","Cases")

    barplot(together, border="white", xlab="Group",
            col = c("#66C2A5","#FC8D62","#8DA0CB"),
            ylab = "Relative Abundances",
            cex.lab = 1.5)

    text(.7,1-1/6,
          "Taxon 1",cex = 1.1)

    text(.7,1-3/6,
         "Taxon 2",cex = 1.1)

    text(.7,1-5/6,
         "Taxon 3",cex = 1.1)

  })

}

# Create Shiny app ----
shinyApp(ui, server)
