library(shiny)

# Sigma and Bandwidth values for Gaussian Kernel Smoothing
sigma_values <- c(1, 2, 5, 8, 10, 15, 20)
bandwidth_values <- seq(0, 0.1, length.out = 11)

# Number of Basis values for Thin Plate Spline Regression
basis_values <- c(20, 40, 60, 80, 100, 120, 140, 160, 180)

# Define UI
ui <- fluidPage(
  titlePanel("HRF Estimation & Inference"),
  
  tabsetPanel(
    id = "session_tabs",
    
    # Gaussian Kernel Smoothing Session
    tabPanel("Gaussian Kernel Smoothing",
             tabsetPanel(
               id = "gks_tabs",
               
               # Results Subsession
               tabPanel("Results",
                        sidebarLayout(
                          sidebarPanel(
                            sliderInput("sigma_gks", "Sigma:", min = 1, max = 20, value = 1, step = 1, ticks = TRUE),
                            sliderInput("bandwidth", "Gaussian Kernel Bandwidth:", min = 0, max = 0.1, value = 0, step = 0.01, ticks = TRUE)
                          ),
                          mainPanel(
                            fluidRow(
                              column(6, imageOutput("gks_image1", height = "250px")),
                              column(6, imageOutput("gks_image2", height = "250px"))
                            ),
                            fluidRow(
                              column(4, imageOutput("gks_image3", height = "250px")),
                              column(4, imageOutput("gks_image4", height = "250px")),
                              column(4, imageOutput("gks_image5", height = "250px"))
                            )
                          )
                        )
               ),
               
               tabPanel("Summary",
                        sidebarLayout(
                          sidebarPanel(
                            sliderInput("sigma_gks_summary", "Sigma:", min = 1, max = 20, value = 1, step = 1, ticks = TRUE)
                          ),
                          mainPanel(
                            fluidRow(
                              column(6, imageOutput("gks_summary_plot1", height = "300px")),
                              column(6, imageOutput("gks_summary_plot2", height = "300px"))
                            )
                          )
                        )
               )
          
             )
    ),
    
    # Thin Plate Spline Regression Session
    tabPanel("Thin Plate Spline Regression",
             tabsetPanel(
               id = "tps_tabs",
               
               # Results Subsession
               tabPanel("Results",
                        sidebarLayout(
                          sidebarPanel(
                            sliderInput("sigma_tps", "Sigma:", min = 1, max = 20, value = 1, step = 1, ticks = TRUE),
                            sliderInput("basis", "Number of Basis:", min = 20, max = 180, value = 20, step = 20, ticks = TRUE)
                          ),
                          mainPanel(
                            fluidRow(
                              column(6, imageOutput("tps_image1", height = "250px")),
                              column(6, imageOutput("tps_image2", height = "250px"))
                            ),
                            fluidRow(
                              column(4, imageOutput("tps_image3", height = "250px")),
                              column(4, imageOutput("tps_image4", height = "250px")),
                              column(4, imageOutput("tps_image5", height = "250px"))
                            )
                          )
                        )
               ),
               
               tabPanel("Summary",
                        sidebarLayout(
                          sidebarPanel(
                            sliderInput("sigma_tps_summary", "Sigma:", min = 1, max = 20, value = 1, step = 1, ticks = TRUE)
                          ),
                          mainPanel(
                            fluidRow(
                              column(6, imageOutput("tps_summary_plot1", height = "300px")),
                              column(6, imageOutput("tps_summary_plot2", height = "300px"))
                            )
                          )
                        )
               )
               
             )
    )
  )
)

# Server Logic
server <- function(input, output, session) {
  # Helper function for Gaussian Kernel Smoothing Results
  generateGKSFilePath <- function(sigma, bandwidth, prefix) {
    closest_sigma <- sigma_values[which.min(abs(sigma_values - sigma))]
    closest_bandwidth <- bandwidth_values[which.min(abs(bandwidth_values - bandwidth))]
    sigma_index <- match(closest_sigma, sigma_values)
    file.path("/Users/user/Documents/JHU/research/HRF_Est/Pain/figure/setting3/GKS", 
              sprintf("%s%d-%.2f.png", prefix, sigma_index, closest_bandwidth))
  }
  
  # Helper function for Thin Plate Spline Regression Results
  generateTPSFilePath <- function(sigma, basis, prefix) {
    closest_sigma <- sigma_values[which.min(abs(sigma_values - sigma))]
    basis_index <- match(basis, basis_values)
    sigma_index <- match(closest_sigma, sigma_values)
    file.path("/Users/user/Documents/JHU/research/HRF_Est/Pain/figure/setting3/TPRS", 
              sprintf("%s%d-%d.png", prefix, sigma_index, basis_index))
  }
  
  # Helper function for Summary File Paths
  generateSummaryFilePath <- function(session_type, sigma, prefix) {
    closest_sigma <- sigma_values[which.min(abs(sigma_values - sigma))]
    sigma_index <- match(closest_sigma, sigma_values)
    
    file.path(
      sprintf("/Users/user/Documents/JHU/research/HRF_Est/Pain/figure/setting3/%s", session_type),
      sprintf("%s%d-.png", prefix, sigma_index)
    )
  }
  
  # Gaussian Kernel Smoothing Summary
  output$gks_summary_plot1 <- renderImage({
    list(
      src = generateSummaryFilePath("GKS", input$sigma_gks_summary, "opt"),
      contentType = "image/png",
      alt = "GKS Summary Plot 1"
    )
  }, deleteFile = FALSE)
  
  output$gks_summary_plot2 <- renderImage({
    list(
      src = generateSummaryFilePath("GKS", input$sigma_gks_summary, "overall_opt"),
      contentType = "image/png",
      alt = "GKS Summary Plot 2"
    )
  }, deleteFile = FALSE)
  
  # Thin Plate Spline Regression Summary
  output$tps_summary_plot1 <- renderImage({
    list(
      src = generateSummaryFilePath("TPRS", input$sigma_tps_summary, "opt"),
      contentType = "image/png",
      alt = "TPS Summary Plot 1"
    )
  }, deleteFile = FALSE)
  
  output$tps_summary_plot2 <- renderImage({
    list(
      src = generateSummaryFilePath("TPRS", input$sigma_tps_summary, "overall_opt"),
      contentType = "image/png",
      alt = "TPS Summary Plot 2"
    )
  }, deleteFile = FALSE)
  
  
  # Gaussian Kernel Smoothing Results
  gks_file_prefixes <- c("betafit", "diff", "pvalue", "FPFN", "cohen")
  for (i in 1:5) {
    local({
      image_index <- i
      prefix <- gks_file_prefixes[image_index]
      output[[paste0("gks_image", image_index)]] <- renderImage({
        list(src = generateGKSFilePath(input$sigma_gks, input$bandwidth, prefix),
             contentType = "image/png", 
             alt = paste("GKS Image", image_index),
             width = 300, height = 250)
      }, deleteFile = FALSE)
    })
  }
  
  
  # Thin Plate Spline Regression Results
  tps_file_prefixes <- c("betafit", "diff", "pvalue", "FPFN", "cohen")
  for (i in 1:5) {
    local({
      image_index <- i
      prefix <- tps_file_prefixes[image_index]
      output[[paste0("tps_image", image_index)]] <- renderImage({
        list(src = generateTPSFilePath(input$sigma_tps, input$basis, prefix),
             contentType = "image/png", 
             alt = paste("TPS Image", image_index),
             width = 300, height = 250)
      }, deleteFile = FALSE)
    })
  }
  
}

# Run the Shiny app
shinyApp(ui = ui, server = server)
