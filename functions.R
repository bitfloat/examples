library(terra)
library(svglite)

.get_distributionType <- function(values) {
  # 0=normal, 1=lognormal, 2=beta, 3=gamma, 4=weibull, 5=poisson, 6=binomial, 7=other

  # check if we have enough data
  if (sum(!is.na(values)) < 3){
    if(all(is.na(values))) return(NA) else return(7) # Not enough data, return "111=other"
  }
  v <- na.omit(values)
  n <- length(v)
  m <- mean(v)
  s <- sd(v)

  if (s == 0) return(0)  # constant values, return "normal"

  # check for discrete/count data
  is_discrete <- all(v == round(v)) && length(unique(v)) < 20

  # check for bounded distributions
  is_bounded_01 <- all(v >= 0 & v <= 1) && min(v) > 0.01 && max(v) < 0.99

  # simple rule-based classification focused on distribution family
  if (is_discrete) {
    if (var(v) > 1.2 * mean(v)) {
      return(6)  # greater variance than poisson suggests binomial
    } else if (all(v >= 0)) {
      return(5)  # count data with variance ≈ mean suggests poisson
    }
  }

  if (is_bounded_01) {
    return(2)  # data in [0,1] range suggests beta
  }

  if (all(v > 0)) {
    # distinguish between positive-only distributions using coefficient of variation
    cv <- s / m

    if (cv < 0.3) {
      return(0)  # low CV - treat as normal
    } else if (cv > 1.0) {
      return(1)  # high CV suggests lognormal
    } else if (m / median(v) > 1.05) {
      return(3)  # right skew with moderate CV suggests gamma
    } else {
      return(4)  # moderate CV suggests weibull
    }
  }

  # default is normal
  return(0)
}

.rebuild_distribution <- function(dist_type, median_val, sd_val,
                                  skew_code = NULL, kurt_code = NULL,
                                  n_points = 1000) {
  # Input validation
  if (is.na(median_val) || is.na(sd_val) || is.na(dist_type)) {
    return(NULL)
  }

  # Prevent negative sd
  sd_val <- max(sd_val, 1e-6)

  # Prevent negative median for distributions that require positive values
  if (dist_type %in% c(1, 2, 3, 4)) {
    median_val <- max(median_val, 1e-6)
  }

  # Define the distribution names for clarity
  dist_names <- c(
    "Normal",       # 0
    "Lognormal",    # 1
    "Beta",         # 2
    "Gamma",        # 3
    "Weibull",      # 4
    "Poisson",      # 5
    "Binomial",     # 6
    "Other"         # 7
  )

  # Convert skew_code to approximate skewness value (midpoint of range)
  skew_values <- c(-2.0, -1.0, -0.3, 0, 0.3, 1.0, 2.75, 5.0)
  skewness <- if(!is.null(skew_code)) skew_values[skew_code + 1] else 0

  # Convert kurt_code to approximate kurtosis value
  kurt_values <- c(-1.5, -0.75, -0.3, 0, 0.3, 1.25, 3.5, 7.0)
  kurtosis <- if(!is.null(kurt_code)) kurt_values[kurt_code + 1] else 0

  # Get distribution name
  dist_name <- dist_names[dist_type + 1]

  # Define x-axis range based on distribution type and kurtosis
  if (dist_type == 0) {
    # Normal: symmetric around median with range adjusted by kurtosis
    x_range <- seq(
      median_val - (4 + 0.5 * max(0, kurtosis)) * sd_val,
      median_val + (4 + 0.5 * max(0, kurtosis)) * sd_val,
      length.out = n_points
    )
  } else if (dist_type %in% c(1, 3, 4)) {
    # Right-skewed distributions: range adjusted by skewness
    left_extend <- max(0, 2 - 0.5 * skewness)
    right_extend <- 4 + skewness
    x_range <- seq(
      max(0, median_val - left_extend * sd_val),
      median_val + right_extend * sd_val,
      length.out = n_points
    )
  } else if (dist_type == 2) {
    # Beta: bounded [0,1]
    x_range <- seq(0, 1, length.out = n_points)
  } else if (dist_type %in% c(5, 6)) {
    # Count distributions: integer range with kurtosis adjustment
    max_val <- ceiling(median_val + (4 + 0.5 * max(0, kurtosis)) * sd_val)
    x_range <- 0:max_val
  } else {
    # Default (Other): similar to normal but wider range
    x_range <- seq(
      median_val - 5 * sd_val,
      median_val + 5 * sd_val,
      length.out = n_points
    )
  }

  # Calculate distribution parameters and density based on distribution type
  if (dist_type == 0) {
    # Normal distribution
    # For normal: median equals mean
    mean_val <- median_val
    density_vals <- dnorm(x_range, mean = mean_val, sd = sd_val)
    params <- list(mean = mean_val, sd = sd_val)

  } else if (dist_type == 1) {
    # Lognormal distribution
    # Use skewness to improve parameter estimation
    if (skewness <= 0.1) {
      # Default method for near-symmetric cases
      meanlog <- log(median_val)
      # CV = sqrt(exp(sdlog^2) - 1)
      sdlog <- sqrt(log(1 + (sd_val/median_val)^2))
    } else {
      # Use skewness for parameter estimation
      # Skewness = (exp(sdlog^2)+2) * sqrt(exp(sdlog^2)-1)
      # Solve for sdlog numerically
      target_func <- function(sdl) {
        (exp(sdl^2)+2) * sqrt(exp(sdl^2)-1) - skewness
      }

      # If target_func(0.1) and target_func(2) have opposite signs, use uniroot
      if (target_func(0.1) * target_func(2) < 0) {
        sdlog <- uniroot(target_func, c(0.1, 2))$root
      } else {
        # Fallback to approximation
        sdlog <- sqrt(log(1 + (skewness/2)^2))
      }

      # With sdlog determined, calculate meanlog to match the median
      meanlog <- log(median_val) - sdlog^2/2
    }

    density_vals <- dlnorm(x_range, meanlog = meanlog, sdlog = sdlog)
    params <- list(meanlog = meanlog, sdlog = sdlog)

  } else if (dist_type == 2) {
    # Beta distribution
    # Use skewness and kurtosis to improve alpha/beta estimation
    if (abs(skewness) < 0.1) {
      # Symmetric case: alpha = beta
      alpha <- 0.5 * ((1/9) * median_val * (1 - median_val) / (sd_val^2) - 1)
      beta <- alpha
    } else {
      # Use skewness for asymmetric cases
      # For Beta: skewness = 2(beta-alpha)*sqrt(alpha+beta+1)/((alpha+beta+2)*sqrt(alpha*beta))
      # Use approximation based on method of moments
      mean_val <- median_val + skewness * sd_val / 3  # Approximation
      mean_val <- max(0.01, min(0.99, mean_val))  # Keep in range

      # Method of moments with mean and variance
      variance <- sd_val^2
      temp <- mean_val * (1 - mean_val) / variance - 1

      if (temp <= 0) {
        # Fallback for invalid parameters
        alpha <- max(0.5, 2 / (1 + abs(skewness)))
        beta <- max(0.5, 2 / (1 + abs(skewness)))
        if (skewness > 0) {
          alpha <- alpha / 2
        } else if (skewness < 0) {
          beta <- beta / 2
        }
      } else {
        alpha <- mean_val * temp
        beta <- (1 - mean_val) * temp
      }
    }

    # Adjust for kurtosis if needed
    if (kurtosis > 1 && alpha > 1 && beta > 1) {
      # Heavy tails: reduce both parameters
      factor <- max(0.5, 1 - kurtosis/10)
      alpha <- alpha * factor
      beta <- beta * factor
    } else if (kurtosis < -0.5 && alpha < 1 && beta < 1) {
      # Light tails: increase both parameters
      factor <- 1 + abs(kurtosis)/3
      alpha <- alpha * factor
      beta <- beta * factor
    }

    # Ensure parameters are at least 0.1
    alpha <- max(0.1, alpha)
    beta <- max(0.1, beta)

    density_vals <- dbeta(x_range, shape1 = alpha, shape2 = beta)
    params <- list(alpha = alpha, beta = beta)

  } else if (dist_type == 3) {
    # Gamma distribution
    # Use skewness for direct parameter estimation
    if (skewness <= 0.1) {
      # Near-symmetric case: default method
      mean_val <- median_val * 1.1
      shape <- (mean_val / sd_val)^2
      scale <- sd_val^2 / mean_val
    } else {
      # For gamma: skewness = 2/sqrt(shape)
      shape <- 4 / (skewness^2)

      # With shape determined, calculate scale to match median
      scale <- median_val / qgamma(0.5, shape = shape, scale = 1)
    }

    # Adjust for kurtosis if needed
    if (!is.null(kurt_code) && kurt_code >= 5) {
      # Heavy-tailed: reduce shape parameter
      shape <- shape * 0.8
      # Recalculate scale to maintain median
      scale <- median_val / qgamma(0.5, shape = shape, scale = 1)
    }

    # Ensure shape parameter is at least 0.1
    shape <- max(0.1, shape)

    density_vals <- dgamma(x_range, shape = shape, scale = scale)
    params <- list(shape = shape, scale = scale)

  } else if (dist_type == 4) {
    # Weibull distribution
    # Use skewness to estimate shape parameter
    if (abs(skewness - 0.6) < 0.3) {
      # Near-standard Weibull
      shape <- 2
    } else if (skewness < 0.4) {
      # More symmetric than standard Weibull
      shape <- 3.6 - skewness
    } else {
      # More skewed than standard Weibull
      shape <- max(0.5, 2 - (skewness - 0.6)/2)
    }

    # With shape determined, calculate scale to match median
    scale <- median_val / (log(2)^(1/shape))

    density_vals <- dweibull(x_range, shape = shape, scale = scale)
    params <- list(shape = shape, scale = scale)

  } else if (dist_type == 5) {
    # Poisson distribution
    # For Poisson: mean = median (approximately)
    lambda <- max(0.1, median_val)

    density_vals <- dpois(x_range, lambda = lambda)
    params <- list(lambda = lambda)

  } else if (dist_type == 6) {
    # Binomial distribution
    # Estimate size and prob
    n_trials <- ceiling(median_val * 4)
    n_trials <- max(10, min(100, n_trials))

    # Calculate p from median and n
    p <- median_val / n_trials
    p <- min(max(0.01, p), 0.99)

    density_vals <- dbinom(x_range, size = n_trials, prob = p)
    params <- list(size = n_trials, prob = p)

  } else {
    # Other/unknown distribution - default to normal
    mean_val <- median_val

    # Adjust standard deviation based on kurtosis
    sd_adjusted <- sd_val
    if (!is.null(kurt_code)) {
      if (kurt_code >= 5) {  # Heavy tails
        sd_adjusted <- sd_val * (1 + 0.2 * (kurt_code - 4))
      } else if (kurt_code <= 2) {  # Light tails
        sd_adjusted <- sd_val * (1 - 0.1 * (3 - kurt_code))
      }
    }

    density_vals <- dnorm(x_range, mean = mean_val, sd = sd_adjusted)
    params <- list(mean = mean_val, sd = sd_adjusted)
  }

  # Return the results
  return(list(
    x = x_range,
    density = density_vals,
    dist_name = dist_name,
    parameters = c(
      list(
        median = median_val,
        sd = sd_val,
        dist_type = dist_type,
        skewness = skewness,
        kurtosis = kurtosis
      ),
      params
    )
  ))
}

.get_skewness <- function(values) {

  # check if we have enough data
  if (sum(!is.na(values)) < 3){
    if(all(is.na(values))) return(NA) else return(7) # Not enough data, return "111=other"
  }

  v <- na.omit(values)
  n <- length(v)
  m <- mean(v)
  s <- sd(v)

  if (s == 0) return(3)  # zero skew for constant values

  # calculate skewness
  skew <- sum((v - m)^3) / (n * s^3)

  # Encode in 3 bits
  if (skew < -1.5) return(0)           # 000=highly negative skewed
  else if (skew < -0.5) return(1)      # 001=moderately negative skewed
  else if (skew < -0.1) return(2)      # 010=slightly negative skewed
  else if (skew < 0.1) return(3)       # 011=approximately symmetric
  else if (skew < 0.5) return(4)       # 100=slightly positive skewed
  else if (skew < 1.5) return(5)       # 101=moderately positive skewed
  else if (skew < 4.0) return(6)       # 110=highly positive skewed
  else return(7)

}

.get_kurtosis <- function(values) {

  # check if we have enough data
  if (sum(!is.na(values)) < 3){
    if(all(is.na(values))) return(NA) else return(7) # Not enough data, return "111=other"
  }
  v <- na.omit(values)
  n <- length(v)
  m <- mean(v)
  s <- sd(v)

  if (s == 0) return(3)  # mesokurtic for constant values

  # calculate excess kurtosis
  kurt <- sum((v - m)^4) / (n * s^4) - 3

  # encode in 3 bits
  if (kurt < -1.0) return(0)           # 000=very light-tailed
  else if (kurt < -0.5) return(1)      # 001=light-tailed
  else if (kurt < -0.1) return(2)      # 010=slightly light-tailed
  else if (kurt < 0.1) return(3)       # 011=mesokurtic (normal-like)
  else if (kurt < 0.5) return(4)       # 100=slightly heavy-tailed
  else if (kurt < 2.0) return(5)       # 101=moderately heavy-tailed
  else if (kurt < 5.0) return(6)       # 110=heavy-tailed
  else return(7)

}

.get_tailThickness <- function(dist_type, kurtosis){

  if (!is.null(kurtosis)) {

    if (kurtosis <= 1) return(0)                 # thin tails (sub-gaussian)
    else if(kurtosis <= 3) return(1)             # normal-like tails
    else if(kurtosis <= 5) return(2)             # moderately heavy tails
    else return(3)                               # heavy tails

  } else {

    # Fall back to distribution type
    if (dist_type %in% c(0, 6)) return(1)        # normal, Binomial: normal-like tails
    else if(dist_type %in% c(1, 3, 4)) return(2) # lognormal, gamma, weibull: moderately heavy tails
    else if(dist_type == 7) return(3)            # other: assume heavy tails to be conservative
    else return(0)                               # beta, others: thin tails

  }

}

.get_tailDirection <- function(dist_type, skewness){

  if (!is.null(skewness)) {

    if (skewness >= 4) return(1)                  # positive skew: right tail risk
    else if(skewness <= 2) return(0)              # negative skew: left tail risk
    else as.integer(dist_type %in% c(1, 3, 4, 5)) # distribution type for symmetric cases

  } else {
    as.integer(dist_type %in% c(1, 3, 4, 5))      # 1 for right-skewed distributions
  }

}

.simulate_uncertaintySource <- function(land, elevation, footprint){

  # assumption that areas with low footprint have worse data cover and thus more uncertainty
  data_coverage <- footprint < 10

  # assumption that topographically more complex areas have more uncertainty
  slope <- terrain(elevation, "slope", unit = "degrees")
  topo_complexity <- focal(slope, w = matrix(1,3,3), fun = var, na.rm = TRUE)
  topo_complexity <- (topo_complexity - min(values(topo_complexity), na.rm = TRUE)) /
    (max(values(topo_complexity), na.rm = TRUE) -
       min(values(topo_complexity), na.rm = TRUE))
  topo_complexity <- as.logical(topo_complexity)

  # assumption that pixels with more landcover classes have more uncertainty
  land_complexity <- 0
  for (layer in names(land)) {
    temp <- focal(land[[layer]] > 0.3, w = matrix(1,3,3), fun = sum, na.rm = TRUE)
    land_complexity <- land_complexity + (temp > 1 & temp < 9)
  }

  # assumption that pixels with different landcover in the neighbourhood have more uncertainty
  edges <- 0
  for (layer in names(land)) {
    temp <- focal(land[[layer]] > 0.3, w = matrix(1,3,3), fun = function(x) {
      if (all(is.na(x))) return(NA)
      any(x, na.rm = TRUE) && !all(x, na.rm = TRUE)
    })
    edges <- edges + temp
  }

  out <-ifel(topo_complexity > 0.3, 1,
             ifel(data_coverage > 0.3, 0,
                   ifel(land_complexity > 0 | edges > 0, 2, 3)))

  names(out) <- "uncertainty_source"
  levels(out) <- data.frame(id = 0:3,
                            uncertainty = c("data coverage", "topography",
                                            "heterogeneity", "edge effects"))

  return(out)
}


.exceedance_probability <- function(animals_mean, animals_sd, dist_values,
                                    skewness, kurtosis, carrying_capacity) {

  # Handle NAs in inputs at the beginning
  carrying_capacity <- ifel(is.na(carrying_capacity), 0, carrying_capacity)
  animals_mean <- ifel(is.na(animals_mean), 0, animals_mean)
  animals_sd <- ifel(is.na(animals_sd) | animals_sd <= 0, 0.01, animals_sd)
  dist_values <- ifel(is.na(dist_values), 0, dist_values)
  skewness <- ifel(is.na(skewness), 3, skewness)
  kurtosis <- ifel(is.na(kurtosis), 3, kurtosis)

  probability <- ifel(animals_mean > carrying_capacity, 1.0, NA)

  # Calculate tail properties for all cells
  tail_thickness <- app(c(dist_values, kurtosis), function(vals) {
    .get_tailThickness(vals[1], vals[2])
  })

  tail_direction <- app(c(dist_values, skewness), function(vals) {
    .get_tailDirection(vals[1], vals[2])
  })

  # Calculate z-scores
  z_scores <- (carrying_capacity - animals_mean) / animals_sd

  # Determine if in risk tail
  in_risk_tail <- ((tail_direction == 1) & (z_scores > 0)) |
    ((tail_direction == 0) & (z_scores < 0))

  # Calculate base probabilities
  abs_z_scores <- abs(z_scores)
  base_probs <- app(abs_z_scores, function(z) pnorm(z, lower.tail = FALSE))

  # Create a raster for the multipliers based on tail thickness
  multipliers <- app(tail_thickness, function(t) {
    c(0.7, 1.0, 1.5, 2.5)[t + 1]
  })

  # Calculate adjusted probabilities for risk tail
  adjusted_probs <- base_probs * multipliers
  adjusted_probs <- ifel(adjusted_probs > 1.0, 1.0, adjusted_probs)

  # Use the appropriate probability based on whether we're in the risk tail
  final_probs <- ifel(in_risk_tail, adjusted_probs, base_probs)
  final_probs <- ifel(animals_mean == 0 & carrying_capacity == 0, 0, final_probs)

  # Update our probability raster where values aren't already set
  probability <- ifel(is.na(probability), final_probs, probability)

  names(probability) <- "exceedance_probability"
  return(probability)

}

ensemble_landscape <- function(ensemble, height_scale = 0.1,
                               color_by = "value", colorscale = "Earth",
                               show_quantiles = TRUE,
                               quantile_levels = c(0.05, 0.95)) {

  # calculate mean layer
  mean_layer <- app(ensemble, mean, na.rm = TRUE)

  # calculate quantile layers if requested
  if (show_quantiles) {
    quantile_layers <- list()
    for (q in quantile_levels) {
      q_layer <- app(ensemble, function(x) quantile(x, q, na.rm = TRUE))
      quantile_layers[[paste0("q", q*100)]] <- q_layer
    }
  }

  # choose coloring statistic
  if (color_by == "variance") {
    color_layer <- app(ensemble, var, na.rm = TRUE)
    color_title <- "Variance"
  } else if (color_by == "cv") {
    sd_layer <- app(ensemble, sd, na.rm = TRUE)
    color_layer <- sd_layer / mean_layer
    color_title <- "Coefficient of Variation"
  } else {
    color_layer <- mean_layer
    color_title <- "Value"
  }

  # extract as matrices
  z_mat <- as.matrix(mean_layer, wide = TRUE)
  color_mat <- as.matrix(color_layer, wide = TRUE)

  # flip matrices if necessary
  if (yres(ensemble) > 0) {
    z_mat <- z_mat[nrow(z_mat):1, ]
    color_mat <- color_mat[nrow(color_mat):1, ]

    # also flip quantile matrices if they exist
    if (show_quantiles) {
      for (q in names(quantile_layers)) {
        quantile_layers[[q]] <- as.matrix(quantile_layers[[q]], wide = TRUE)[nrow(as.matrix(quantile_layers[[q]], wide = TRUE)):1, ]
      }
    }
  }

  # Get the extent of the raster
  r_extent <- ext(ensemble)

  # create x and y coordinate vectors (longitude and latitude)
  x_res <- xres(ensemble)
  y_res <- yres(ensemble)

  # get coordinate vectors for cell centers
  x_coords <- seq(from = r_extent[1] + x_res/2, to = r_extent[2] - x_res/2, by = x_res)
  y_coords <- seq(from = r_extent[3] + y_res/2, to = r_extent[4] - y_res/2, by = y_res)

  # make sure the coordinates match the matrix dimensions
  if (length(x_coords) != ncol(z_mat)) {
    x_coords <- seq(from = r_extent[1], to = r_extent[2], length.out = ncol(z_mat))
  }
  if (length(y_coords) != nrow(z_mat)) {
    y_coords <- seq(from = r_extent[3], to = r_extent[4], length.out = nrow(z_mat))
  }

  # calculate aspect ratio for preserving original dimensions
  x_range <- max(x_coords) - min(x_coords)
  y_range <- max(y_coords) - min(y_coords)
  z_range <- max(z_mat, na.rm = TRUE) - min(z_mat, na.rm = TRUE)

  # normalize to make the largest dimension 1.0
  max_range <- max(x_range, y_range)
  x_aspect <- x_range / max_range
  y_aspect <- y_range / max_range
  z_aspect <- (z_range / max_range) * height_scale

  # create the 3D plot
  p <- plot_ly()

  # add quantile surfaces first (so they appear beneath the mean)
  if (show_quantiles) {
    # sort quantiles to process them from most extreme to least extreme
    sorted_quantiles <- sort(quantile_levels)

    # determine base colors based on colorscale
    if (colorscale == "Earth") {
      base_color_lower <- "rgba(36, 67, 15, {alpha})"   # Earth dark green
      base_color_upper <- "rgba(209, 154, 102, {alpha})" # Earth light brown
    } else if (colorscale == "Viridis") {
      base_color_lower <- "rgba(68, 1, 84, {alpha})"   # Viridis dark purple
      base_color_upper <- "rgba(253, 231, 37, {alpha})" # Viridis yellow
    } else if (colorscale == "Plasma") {
      base_color_lower <- "rgba(13, 8, 135, {alpha})"   # Plasma dark blue
      base_color_upper <- "rgba(240, 249, 33, {alpha})" # Plasma yellow
    } else if (colorscale == "Inferno") {
      base_color_lower <- "rgba(0, 0, 4, {alpha})"     # Inferno black
      base_color_upper <- "rgba(252, 255, 164, {alpha})" # Inferno light yellow
    } else if (colorscale == "Cividis") {
      base_color_lower <- "rgba(0, 32, 76, {alpha})"    # Cividis dark blue
      base_color_upper <- "rgba(255, 236, 156, {alpha})" # Cividis light yellow
    } else if (colorscale == "Turbo") {
      base_color_lower <- "rgba(48, 18, 59, {alpha})"   # Turbo dark purple
      base_color_upper <- "rgba(248, 30, 7, {alpha})"   # Turbo red
    } else {
      # default fallback colors
      base_color_lower <- "rgba(0, 0, 255, {alpha})"    # Blue
      base_color_upper <- "rgba(255, 0, 0, {alpha})"    # Red
    }

    # create color list with varying alpha (saturation proxy)
    q_colors <- list()
    midpoint <- 0.5

    for (q in sorted_quantiles) {
      # determine if it's lower or upper quantile
      is_lower <- q < midpoint

      # calculate alpha based on distance from median (0.2-0.5 range)
      alpha <- 0.2 + (abs(q - midpoint) * 0.6)

      # format color string with proper alpha
      if (is_lower) {
        color <- gsub("\\{alpha\\}", alpha, base_color_lower)
      } else {
        color <- gsub("\\{alpha\\}", alpha, base_color_upper)
      }

      q_colors[[paste0("q", q*100)]] <- color
    }

    # add each quantile surface
    for (q in sorted_quantiles) {
      q_name <- paste0("q", q*100)
      q_mat <- as.matrix(quantile_layers[[q_name]], wide = TRUE)

      # get appropriate color
      q_color <- q_colors[[q_name]]

      # add surface
      p <- p %>% add_surface(
        x = x_coords,
        y = y_coords,
        z = q_mat,
        colorscale = list(c(0, 1), c(q_color, q_color)),
        showscale = FALSE,
        opacity = 1,  # Opacity is already in the RGBA colors
        hoverinfo = "text",
        text = paste0(q*100, "% quantile"),
        name = paste0(q*100, "%")
      )
    }
  }

  # add main surface on top
  p <- p %>% add_surface(
    x = x_coords,
    y = y_coords,
    z = z_mat,
    surfacecolor = color_mat,
    colorscale = colorscale,
    colorbar = list(title = color_title),
    name = "Mean"
  )

  # define top-down camera view
  top_down_camera <- list(
    eye = list(x = 0, y = 0, z = 2.5),  # Camera position
    center = list(x = 0, y = 0, z = 0),  # Point camera looks at
    up = list(x = 0, y = 1, z = 0)  # Up direction
  )

  # finalize layout with top-down view and preserved aspect ratio
  p <- p %>% layout(
    scene = list(
      aspectmode = "manual",
      aspectratio = list(x = y_aspect, y = x_aspect, z = z_aspect),
      xaxis = list(title = "Longitude"),
      yaxis = list(title = "Latitude"),
      zaxis = list(title = color_title),
      camera = top_down_camera
    ),
    title = "Ensemble Data Surface",
    legend = list(x = 0, y = 1)  # Move legend to top-left
  )

  # add configuration to move modebar to the left
  p <- p %>% config(
    modeBarButtonsToRemove = c("sendDataToCloud", "editInChartStudio", "lasso2d", "select2d"),
    displaylogo = FALSE,
    toImageButtonOptions = list(
      format = "png",
      filename = "ensemble_plot",
      width = 1200,
      height = 800
    )
  )

  return(p)
}

.export_svg <- function(raster){

  svglite(paste0(getwd(), "/figures/", names(raster), ".svg"),
          width = 10.2, height = 7.8,
          pointsize = 12, bg = "transparent")

  plot(raster)

  dev.off()

}
