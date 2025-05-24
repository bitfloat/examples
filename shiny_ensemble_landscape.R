library(shiny)
library(plotly)
library(ggplot2)
library(dplyr)
library(tidyr)
library(terra)
library(viridis)

# check ensemble object is available
if(is.null(tryCatch(get("ensemble"), error = function(e) NULL))) stop("ensemble stack is missing")

# define colours
myCols <- list("#21908C", "#fde725", "#440154", "white", "black", "lightblue")
names(myCols) <- c("mean", "upper", "lower", "white", "black", "dist")


# Notes for Shiny App Enhancements:
#
#   1. Manual Coordinate Input
#
# - Add text input fields for longitude and latitude
# - Add a "Go to Location" button that updates the current selection
# - Implement a reactive value to track whether selection comes from click or input fields
#
#   2. Custom Transects
#
# - Add functionality to define start and end points
# - Implement transect calculation along arbitrary lines
# - Create new plot to show data along transect

# ui.R
ui <- fluidPage(
  titlePanel("Ensemble Data Explorer"),

  # Top row: 3D plot (larger) and density plot (smaller)
  fluidRow(
    column(8, plotlyOutput("surface3d", height = "600px")),
    column(2, plotOutput("densityPlot", height = "600px"))
  ),

  # Bottom row: X and Y cross-sections side by side
  fluidRow(
    uiOutput("crossSectionRow")
  )
)

# server.R
server <- function(input, output, session) {
  # Calculate raster dimensions and proportions
  raster_dimensions <- reactive({
    n_x <- ncol(ensemble)
    n_y <- nrow(ensemble)

    # Calculate relative width proportions
    total_width <- 12  # Bootstrap total columns

    # Use the ratio of dimensions to determine column widths
    ratio <- n_x / n_y

    # Calculate column widths based on ratio
    # Make sure they add up to total_width
    if (ratio >= 1) {
      # X dimension is larger
      x_cols <- round(total_width * (ratio / (1 + ratio)))
      y_cols <- total_width - x_cols
    } else {
      # Y dimension is larger
      y_cols <- round(total_width * (1 / (1 + ratio)))
      x_cols <- total_width - y_cols
    }

    # Make sure neither gets 0 columns
    x_cols <- max(x_cols, 3)
    y_cols <- max(y_cols, 3)

    # Adjust to fit within 12 columns
    if (x_cols + y_cols > 12) {
      scale_factor <- 12 / (x_cols + y_cols)
      x_cols <- floor(x_cols * scale_factor)
      y_cols <- 12 - x_cols
    }

    return(list(
      n_x = n_x,
      n_y = n_y,
      x_cols = x_cols,
      y_cols = y_cols
    ))
  })

  # Dynamically generate cross-section row UI based on raster dimensions
  output$crossSectionRow <- renderUI({
    dims <- raster_dimensions()

    # Fill the entire available height (using 100% of container)
    fluidRow(
      column(dims$x_cols, div(style = "height: 200px; margin-bottom: 0px; padding: 0px;",
                              plotOutput("xCrossSection", height = "100%"))),
      column(dims$y_cols, div(style = "height: 200px; margin-bottom: 0px; padding: 0px;",
                              plotOutput("yCrossSection", height = "100%")))
    )
  })

  # Quantile levels to display
  quantile_levels <- c(0.1, 0.25, 0.75, 0.9)

  # Store spatial coordinates for cross-section conversion
  spatial_coords <- reactive({
    # Get the extent of the raster
    r_extent <- ext(ensemble)

    # Create x and y coordinate vectors
    x_res <- xres(ensemble)
    y_res <- yres(ensemble)

    # Get coordinate vectors for cell centers
    x_coords <- seq(from = r_extent[1] + x_res/2, to = r_extent[2] - x_res/2, by = x_res)
    y_coords <- seq(from = r_extent[3] + y_res/2, to = r_extent[4] - y_res/2, by = y_res)

    # Make sure the coordinates match the dimensions
    mean_mat <- as.matrix(app(ensemble, mean, na.rm = TRUE), wide = TRUE)
    if (length(x_coords) != ncol(mean_mat)) {
      x_coords <- seq(from = r_extent[1], to = r_extent[2], length.out = ncol(mean_mat))
    }
    if (length(y_coords) != nrow(mean_mat)) {
      y_coords <- seq(from = r_extent[3], to = r_extent[4], length.out = nrow(mean_mat))
    }

    # Return coordinates
    list(x = x_coords, y = y_coords)
  })

  # Data matrices and values (reactive to ensure they're only computed once)
  all_data <- reactive({
    # Mean matrix
    mean_mat <- as.matrix(app(ensemble, mean, na.rm = TRUE), wide = TRUE)

    # Quantile matrices
    q_matrices <- list()
    for (q in quantile_levels) {
      q_layer <- app(ensemble, function(x) quantile(x, q, na.rm = TRUE))
      q_matrices[[paste0("q", q*100)]] <- as.matrix(q_layer, wide = TRUE)
    }

    # Extract all values for density plot
    all_values <- as.vector(values(ensemble))
    all_values <- all_values[!is.na(all_values)]

    # Flip matrices if necessary
    if (yres(ensemble) > 0) {
      mean_mat <- mean_mat[nrow(mean_mat):1, ]
      for (q in names(q_matrices)) {
        q_matrices[[q]] <- q_matrices[[q]][nrow(q_matrices[[q]]):1, ]
      }
    }

    # Return data in a list
    return(list(
      mean = mean_mat,
      quantiles = q_matrices,
      all_values = all_values
    ))
  })

  # Debounce functionality for hover events
  debounced_pos <- reactive({
    hover_data <- event_data("plotly_hover", source = "surface3d")

    if (is.null(hover_data) || is.null(hover_data$x) || is.null(hover_data$y)) {
      return(NULL)
    }

    # Find closest coordinate indices
    coords <- spatial_coords()
    x_idx <- which.min(abs(coords$x - hover_data$x[1]))
    y_idx <- which.min(abs(coords$y - hover_data$y[1]))

    list(
      x = x_idx,
      y = y_idx,
      lon = hover_data$x[1],
      lat = hover_data$y[1]
    )
  }) %>% debounce(500)  # 500ms debounce

  # Store cross-section coordinates
  crossSection <- reactiveValues(x = NULL, y = NULL, lon = NULL, lat = NULL)

  # Update cross-section values when debounced position changes
  observe({
    pos <- debounced_pos()
    if (!is.null(pos)) {
      crossSection$x <- pos$x
      crossSection$y <- pos$y
      crossSection$lon <- pos$lon
      crossSection$lat <- pos$lat
    }
  })

  # Extract cell-specific values from raster ensemble
  cell_distribution <- reactive({
    req(crossSection$x, crossSection$y)

    # Extract cell coordinates
    coords <- spatial_coords()
    x_coord <- coords$x[crossSection$x]
    y_coord <- coords$y[crossSection$y]

    # Extract cell values from all layers of the ensemble
    cell_values <- terra::extract(ensemble, matrix(c(x_coord, y_coord), ncol=2))

    # Remove NA values
    if (is.data.frame(cell_values)) {
      cell_values <- unlist(cell_values)
    }
    cell_values <- cell_values[!is.na(cell_values)]

    return(cell_values)
  })

  # Main 3D plot
  output$surface3d <- renderPlotly({
    # Create enhanced surface plot with spatial coordinates and left controls
    p <- ensemble_landscape(
      ensemble,
      quantile_levels = quantile_levels,
    )

    # Set source for hover events
    p$x$source <- "surface3d"

    p
  })

  # Density plot in top right - now showing cell-specific distribution
  output$densityPlot <- renderPlot({
    # Only render if we have valid coordinates
    req(crossSection$x, crossSection$y)
    req(all_data())

    data <- all_data()
    mean_mat <- data$mean
    q_mats <- data$quantiles

    # Get cell-specific distribution
    cell_values <- cell_distribution()

    # Validate that indices are within bounds
    valid_x <- min(max(1, crossSection$x), ncol(mean_mat))
    valid_y <- min(max(1, crossSection$y), nrow(mean_mat))

    # Get current point statistics
    current_point <- list(
      mean = mean_mat[valid_y, valid_x],
      q10 = q_mats$q10[valid_y, valid_x],
      q25 = q_mats$q25[valid_y, valid_x],
      q75 = q_mats$q75[valid_y, valid_x],
      q90 = q_mats$q90[valid_y, valid_x]
    )

    # Check if we have enough data points for density estimation
    use_density <- length(cell_values) >= 10

    # Create the plot with a more compact design
    p <- ggplot()

    if (use_density) {
      # Use density plot if enough data points
      p <- p + geom_density(aes(x = cell_values),
                            fill = myCols$dist, alpha = 0.6)
    } else {
      # Use histogram if few data points
      p <- p + geom_histogram(aes(x = cell_values),
                              fill = myCols$dist, alpha = 0.6,
                              bins = max(5, length(cell_values)/2))
    }

    # Add vertical lines for the quantiles
    p <- p +
      geom_vline(xintercept = current_point$q10,
                 color = myCols$lower, linetype = "dashed", alpha = 0.7) +
      geom_vline(xintercept = current_point$q25,
                 color = myCols$lower, linetype = "dashed", alpha = 0.9) +
      geom_vline(xintercept = current_point$q75,
                 color = myCols$upper, linetype = "dashed", alpha = 0.9) +
      geom_vline(xintercept = current_point$q90,
                 color = myCols$upper, linetype = "dashed", alpha = 0.7) +
      geom_vline(xintercept = current_point$mean,
                 color = myCols$mean, size = 1.2)

    # Add annotations (more compact)
    p <- p +
      annotate("text", x = current_point$mean, y = Inf,
               label = "μ", color = myCols$mean, hjust = -0.2, vjust = 2, size = 3)

    # Add sample size info
    p <- p +
      annotate("text", x = -Inf, y = -Inf,
               label = paste("n =", length(cell_values)),
               hjust = -0.1, vjust = -0.5, size = 2.5)

    # Customize theme for a more compact display
    p <- p +
      theme_minimal() +
      labs(
        title = "Cell Distribution",
        subtitle = paste0("(", round(crossSection$lon, 3), ", ",
                          round(crossSection$lat, 3), ")"),
        x = "Value",
        y = if(use_density) "Density" else "Count"
      ) +
      theme(
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.margin = margin(5, 5, 5, 5)
      )

    p
  })

  # X cross-section plot (showing values along a constant y)
  output$xCrossSection <- renderPlot({
    # Only render if we have valid y coordinate
    req(crossSection$y)
    req(all_data())
    req(spatial_coords())

    data <- all_data()
    mean_mat <- data$mean
    q_mats <- data$quantiles
    coords <- spatial_coords()

    # Validate that the index is within bounds
    valid_y <- min(max(1, crossSection$y), nrow(mean_mat))

    # Use longitude values for x-axis
    x_values <- coords$x
    mean_data <- mean_mat[valid_y, ]

    # Prepare data for ggplot
    plot_data <- data.frame(
      x = x_values,
      mean = mean_data
    )

    # Add quantile data
    for (q in quantile_levels) {
      q_name <- paste0("q", q*100)
      plot_data[[q_name]] <- q_mats[[q_name]][valid_y, ]
    }

    # Create the ggplot
    p <- ggplot(plot_data, aes(x = x)) +
      # Add ribbon for 10-90% interval
      geom_ribbon(aes(ymin = q10, ymax = q90),
                  fill = myCols$lower, alpha = 0.15) +

      # Add ribbon for 25-75% interval (interquartile range)
      geom_ribbon(aes(ymin = q25, ymax = q75),
                  fill = myCols$lower, alpha = 0.25) +

      # Add quantile lines
      geom_line(aes(y = q10, color = "10%"),
                linetype = "dashed", size = 0.5) +
      geom_line(aes(y = q25, color = "25%"),
                linetype = "dashed", size = 0.75) +
      geom_line(aes(y = q75, color = "75%"),
                linetype = "dashed", size = 0.75) +
      geom_line(aes(y = q90, color = "90%"),
                linetype = "dashed", size = 0.5) +

      # Add mean line on top
      geom_line(aes(y = mean, color = "Mean"),
                size = 1.2) +

      # If we have a current longitude position, highlight it
      geom_vline(xintercept = crossSection$lon,
                 color = myCols$black, linetype = "solid", alpha = 0.5) +

      # Custom colors
      scale_color_manual(
        name = "Statistics",
        values = c("Mean" = myCols$mean,
                   "10%" = adjustcolor(myCols$lower, alpha.f = 0.8),
                   "25%" = adjustcolor(myCols$lower, alpha.f = 0.9),
                   "75%" = adjustcolor(myCols$upper, alpha.f = 0.9),
                   "90%" = adjustcolor(myCols$upper, alpha.f = 0.8))
      ) +

      # Customize theme
      theme_minimal() +
      labs(
        title = paste("Longitude Cross-section at latitude =", round(crossSection$lat, 3)),
        x = "Longitude",
        y = "Value"
      ) +
      theme(
        legend.position = "top",
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        legend.key.size = unit(0.5, "lines"),  # Smaller legend
        plot.margin = unit(c(0, 0, 0, 0), "mm")  # Minimize margins
      )

    p
  })

  # Y cross-section plot (showing values along a constant x)
  output$yCrossSection <- renderPlot({
    # Only render if we have valid x coordinate
    req(crossSection$x)
    req(all_data())
    req(spatial_coords())

    data <- all_data()
    mean_mat <- data$mean
    q_mats <- data$quantiles
    coords <- spatial_coords()

    # Validate that the index is within bounds
    valid_x <- min(max(1, crossSection$x), ncol(mean_mat))

    # Use latitude values for y-axis
    y_values <- coords$y
    mean_data <- mean_mat[, valid_x]

    # Prepare data for ggplot
    plot_data <- data.frame(
      y = y_values,
      mean = mean_data
    )

    # Add quantile data
    for (q in quantile_levels) {
      q_name <- paste0("q", q*100)
      plot_data[[q_name]] <- q_mats[[q_name]][, valid_x]
    }

    # Create the ggplot
    p <- ggplot(plot_data, aes(x = y)) +
      # Add ribbon for 10-90% interval
      geom_ribbon(aes(ymin = q10, ymax = q90),
                  fill = myCols$lower, alpha = 0.15) +

      # Add ribbon for 25-75% interval (interquartile range)
      geom_ribbon(aes(ymin = q25, ymax = q75),
                  fill = myCols$lower, alpha = 0.25) +

      # Add quantile lines
      geom_line(aes(y = q10, color = "10%"),
                linetype = "dashed", size = 0.5) +
      geom_line(aes(y = q25, color = "25%"),
                linetype = "dashed", size = 0.75) +
      geom_line(aes(y = q75, color = "75%"),
                linetype = "dashed", size = 0.75) +
      geom_line(aes(y = q90, color = "90%"),
                linetype = "dashed", size = 0.5) +

      # Add mean line on top
      geom_line(aes(y = mean, color = "Mean"),
                size = 1.2) +

      # If we have a current latitude position, highlight it
      geom_vline(xintercept = crossSection$lat,
                 color = myCols$black, linetype = "solid", alpha = 0.5) +

      # Custom colors
      scale_color_manual(
        name = "Statistics",
        values = c("Mean" = myCols$mean,
                   "10%" = adjustcolor(myCols$lower, alpha.f = 0.8),
                   "25%" = adjustcolor(myCols$lower, alpha.f = 0.9),
                   "75%" = adjustcolor(myCols$upper, alpha.f = 0.9),
                   "90%" = adjustcolor(myCols$upper, alpha.f = 0.8))
      ) +

      # Customize theme
      theme_minimal() +
      labs(
        title = paste("Latitude Cross-section at longitude =", round(crossSection$lon, 3)),
        x = "Latitude",
        y = "Value"
      ) +
      theme(
        legend.position = "none",  # Hide legend on second plot
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 9),
        axis.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), "mm")  # Minimize margins
      )

    p
  })
}

shinyApp(ui, server)
