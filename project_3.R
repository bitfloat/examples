library(terra)
library(geodata)
library(tidyverse)
library(bitfield)

set.seed(42)

xmin <- 7.6
xmax <- 9.3
ymin <- 47.5
ymax <- 48.8
area_extent <- ext(xmin, xmax, ymin, ymax)


# 1. get (bio)geophysical basis ----
#
## get footprint
footprint <- footprint(year = 2009, path = tempdir()) |>
  crop(area_extent)
names(footprint) <- "human_footprint"


# 2. decode bitfield from project_1 to retrieve layers ----
# animals_mean
# animals_sd
# landcover


# 3. decode bitfield from project_2 to retrieve layers ----
# carrying_capacity
# exceedance
# exceedance_risk
# deficit_type
# deficit_magnitude


# 4. simulate market distortions ----
## Distortion type ----
# 0=none, 1=subsidy, 2=tax, 3=regulation, 4=externality, 5=monopoly, 6=information, 7=public good
distortion_type <- footprint

# determine landscape context
is_agriculture <- landcover$cropland > 0.3
is_forest <- landcover$trees > 0.3
is_grassland <- landcover$grassland > 0.3
is_urban_proximate <- footprint > 15
is_protected_area <- footprint < 8 & (landcover$trees > 0.5 | landcover$water > 0.2)
is_remote <- footprint < 12 & !is_urban_proximate

# assign values based on heuristics
values(distortion_type) <- 0                                  # no distortion (base state)
distortion_type[is_agriculture & !is_protected_area] <- 1     # subsidy, common in productive agricultural areas
distortion_type[is_protected_area & is_urban_proximate] <- 2  # tax, sensitive transition zones
distortion_type[is_forest & is_urban_proximate] <- 3          # regulation, natural areas with moderate human activity
distortion_type[is_agriculture & is_urban_proximate &
                  !is_protected_area] <- 4                    # externality, intensive use areas without compensating mechanisms
distortion_type[is_remote & is_agriculture] <- 5              # monopoly, isolated areas with concentrated economic activity
distortion_type[is_grassland & is_forest] <- 6                # information failure, transition areas between different land use types
distortion_type[is_protected_area & !is_urban_proximate] <- 7 # public good, providing ecosystem services

names(distortion_type) <- "distortion_type"
levels(distortion_type) <- data.frame(id = 0:7,
                                      distortion = c("none", "subsidy", "tax",
                                                     "regulation", "externality",
                                                     "monopoly", "information",
                                                     "public good"))

## distortion magnitude ----
distortion_magnitude <- footprint
values(distortion_magnitude) <- 0

# normalize to 0-1 scale
footprint_norm <- (footprint - min(values(footprint), na.rm=TRUE)) /
  (max(values(footprint), na.rm=TRUE) - min(values(footprint), na.rm=TRUE))

distortion_magnitude <- ifel(distortion_type == 1,
                             round(2 + 2 * footprint_norm),
                             distortion_magnitude) # subsidy; at least distortion of 2 for subsidies and scaling factor of 2 as well
distortion_magnitude <- ifel(distortion_type == 2,
                             round(1 + 2 * footprint_norm),
                             distortion_magnitude) # tax
distortion_magnitude <- ifel(distortion_type == 3,
                             round(1 + 3 * footprint_norm),
                             distortion_magnitude) # regulation
distortion_magnitude <- ifel(distortion_type == 4,
                             round(3 + 3 * footprint_norm),
                             distortion_magnitude) # externality
distortion_magnitude <- ifel(distortion_type == 5,
                             round(4 + 3 * footprint_norm),
                             distortion_magnitude) # monopoly
distortion_magnitude <- ifel(distortion_type == 6,
                             round(2 + 3 * footprint_norm),
                             distortion_magnitude) # information
distortion_magnitude <- ifel(distortion_type == 7,
                             round(2 + 4 * footprint_norm),
                             distortion_magnitude) # public good; very high scaling as undersupply worsens dramatically with intensity

names(distortion_magnitude) <- "distortion_magnitude"
levels(distortion_magnitude) <- data.frame(id = 0:7,
                                           magnitude = c("0-5%", "5-10%", "10-20%",
                                                          "20-40%", "40-70%",
                                                          "70-100%", ">100%",
                                                          "variable"))

# plot
plot(c(lc_dom, footprint, distortion_type, distortion_magnitude))


# 5. calculate GDP adjustment estimates ----
# estimate hidden ecological costs as % of conventional GDP
exceedance_norm <- exceedance / 10  # Normalize to 0-1 scale

# water and soil limitations typically have higher remediation costs
resource_weights <- classify(deficit_type,
                             rcl = matrix(c(0, 0.5,    # no limitation: moderate impact
                                            1, 1.3,    # water limitation: high cost
                                            2, 0.8,    # forage limitation: moderate cost
                                            3, 1.2),   # soil limitation: high cost
                                          ncol=2, byrow=TRUE))

# resource adjustment component
resource_adjustment <- resource_weights * (deficit_magnitude / 7)

# market distortion component
distortion_weights <- classify(distortion_type,
                               rcl = matrix(c(0, 0.2,    # no distortion: minimal impact
                                              1, 1.1,    # subsidy: often hides true costs
                                              2, 0.8,    # tax: may already internalize costs
                                              3, 0.9,    # regulation: partial internalization
                                              4, 1.3,    # externality: significant hidden costs
                                              5, 1.2,    # monopoly: market inefficiency costs
                                              6, 1.0,    # information failure: moderate impact
                                              7, 0.7),   # public good: often undervalued
                                            ncol=2, byrow=TRUE))

# distortion effect
distortion_effect <- distortion_weights * (distortion_magnitude / 7)

# calculate combined GDP adjustment percentage with weighted components
gdp_adjustment <-
  (0.4 * exceedance_norm) +           # highest weight: carrying capacity exceedance
  (0.3 * resource_adjustment) +       # medium weight: resource limitations
  (0.3 * distortion_effect)           # medium weight: market distortions

# ecological systems often show non-linear impacts beyond certain thresholds
gdp_adjustment <- ifel(exceedance > 3,
                       gdp_adjustment * 1.5,  # 50% amplification for severe exceedance
                       gdp_adjustment)

# Constrain values between 0-1 (0-100%)
gdp_adjustment <- clamp(gdp_adjustment, 0, 1)

# Scale to 0-63 range for 6-bit encoding
gdp_adjustment <- round(gdp_adjustment * 63)
names(gdp_adjustment) <- "gdp_adjustment"


# 6. calculate economic-ecological misalignment index ----
# component 1: ecological pressure - how much current use exceeds sustainable capacity
ecological_pressure <- exceedance / 10  # normalize to 0-1 scale, capped at 1

# component 2: resource constraint severity - limitations in critical resources
resource_constraint <- deficit_magnitude / 7

# component 3: market failure - economic signals that fail to reflect ecological realities
market_failure <- distortion_magnitude / 7

# component 4: hidden costs - ecological costs not captured in economic accounting
hidden_costs <- gdp_adjustment / 63

# base misalignment components
base_misalignment <-
  (0.35 * ecological_pressure) +    # highest weight: direct ecological overshoot
  (0.25 * resource_constraint) +    # medium weight: resource limitations
  (0.25 * market_failure) +         # medium weight: market distortions
  (0.15 * hidden_costs)             # lower weight: derived metric

# calculate neighborhood effect to capture spatial interdependencies that
# reflect how sustainability issues cross property boundaries and require
# landscape-level solutions
neighborhood_effect <- focal(base_misalignment,
                             w = matrix(1, 3, 3),  # 3x3 moving window
                             fun = "mean",
                             na.rm = TRUE)

# normalize neighborhood effect
ne_min <- global(neighborhood_effect, "min", na.rm = TRUE)[1, 1]
ne_max <- global(neighborhood_effect, "max", na.rm = TRUE)[1, 1]
neighborhood_norm <- (neighborhood_effect - ne_min) / (ne_max - ne_min)

# calculate final misalignment with spatial effects
misalignment <- (0.7 * base_misalignment) + (0.3 * neighborhood_norm)

# constrain to 0-1 range
misalignment <- clamp(misalignment, 0, 1)

# scale to 0-255 range for 8-bit encoding
misalignment <- round(misalignment * 255)
names(misalignment) <- "economic-ecological_misalignment"


# 7. generate system trajectory classifications ----
# this adds a critical temporal dimension to the analysis

# normalize inputs for clarity
misalign_norm <- misalignment/255
deficit_norm <- deficit_magnitude/7

# generate individual condition masks to improve readability
is_sustainable <- exceedance < 1
is_moderate_exceedance <- exceedance >= 1 & exceedance < 1.5
is_high_exceedance <- exceedance >= 1.5
is_approaching_threshold <- exceedance >= 0.9 & exceedance < 1

is_low_misalignment <- misalign_norm < 0.3
is_moderate_misalignment <- misalign_norm >= 0.3 & misalign_norm < 0.6
is_high_misalignment <- misalign_norm >= 0.6

is_resource_limited <- deficit_norm > 0.5

system_trajectory <- carrying_capacity
values(system_trajectory) <- NA

# stable systems: low misalignment, well below carrying capacity
system_trajectory <- ifel(misalign_norm < 0.3 & exceedance < 0.7, 0, system_trajectory)

# sustainable intensification: low misalignment, approaching but not exceeding capacity
system_trajectory <- ifel(misalign_norm < 0.3 & exceedance >= 0.7 & exceedance < 1, 1, system_trajectory)

# approaching capacity: moderate misalignment, close to capacity
system_trajectory <- ifel(misalign_norm >= 0.3 & misalign_norm < 0.4 & exceedance >= 0.8 & exceedance < 1, 2, system_trajectory)

# cyclic fluctuation: significant resource limitations driving oscillations
system_trajectory <- ifel(deficit_magnitude/7 > 0.6 & exceedance >= 0.7 & exceedance <= 1.3, 3, system_trajectory)

# gradual degradation: moderate misalignment, moderately exceeding capacity
system_trajectory <- ifel(misalign_norm >= 0.3 & misalign_norm < 0.6 & exceedance >= 1 & exceedance < 1.5, 4, system_trajectory)

# rapid degradation: high misalignment, significantly exceeding capacity
system_trajectory <- ifel(misalign_norm >= 0.6 & exceedance >= 1.5, 5, system_trajectory)

# approaching threshold: moderate to high misalignment, near critical threshold
system_trajectory <- ifel(misalign_norm >= 0.4 & exceedance >= 0.9 & exceedance < 1.1, 6, system_trajectory)

# post-threshold reorganization: very high misalignment but reduced pressure
system_trajectory <- ifel(misalign_norm > 0.7 & exceedance < 0.9, 7, system_trajectory)

# default to gradual degradation if exceeding capacity
system_trajectory <- ifel(is.na(system_trajectory) & exceedance >= 1, 4, system_trajectory)

# default to stable for any remaining areas
system_trajectory <- ifel(is.na(system_trajectory), 0, system_trajectory)

names(system_trajectory) <- "system_trajectory"
levels(system_trajectory) <- data.frame(
  id = 0:7,
  trajectory = c("stable", "sustainable intensification", "approaching capacity",
    "cyclic fluctuation", "gradual degradation", "rapid degradation",
    "approaching threshold", "post-threshold reorganization")
)

# plot
plot(c(gdp_adjustment, misalignment, system_trajectory))


# 9. calculate intervention priorities ----
#
# different trajectories have different implications for intervention urgency
trajectory_weights <- classify(system_trajectory,
                               rcl = matrix(c(
                                 0, 0.0,  # stable: lowest priority
                                 1, 0.2,  # sustainable intensification: low priority
                                 2, 0.3,  # approaching capacity: moderate priority
                                 3, 0.5,  # cyclic fluctuation: moderate-high priority
                                 4, 0.7,  # gradual degradation: high priority
                                 5, 0.9,  # rapid degradation: very high priority
                                 6, 1.0,  # approaching threshold: highest priority
                                 7, 0.8   # post-threshold reorganization: high priority (window of opportunity)
                               ), ncol=2, byrow=TRUE))

# higher values indicate areas where intervention is more critical
intervention_priority <-
  (0.35 * (misalignment/255)) +                      # economic-ecological misalignment (35%)
  (0.30 * min(exceedance/10, 1)) +                   # carrying capacity exceedance (30%), capped at 1
  (0.15 * (deficit_magnitude/7)) +                   # resource limitation severity (15%)
  (0.20 * trajectory_weights)                        # system trajectory (increased to 20%)

# Scale to 0-31 range for 5-bit encoding
intervention_priority <- round(intervention_priority * 31)
names(intervention_priority) <- "intervention_priority"


# 10. calculate cross-sectoral synergy potential ----
# identifies opportunities where food security interventions
# simultaneously advance other development goals
water_synergy <- deficit_type == 1
ecosystem_synergy <- exceedance > 1.5 & system_trajectory >= 4
economic_synergy <- misalignment/255 > 0.6
climate_synergy <- deficit_magnitude/7 > 0.5

# combine into a single bitmap (0-15)
cross_sectoral_synergy <-
  (water_synergy * 8) +          # bit 3: water security
  (ecosystem_synergy * 4) +      # bit 2: ecosystem conservation
  (economic_synergy * 2) +       # bit 1: economic development
  (climate_synergy * 1)          # bit 0: climate resilience
names(cross_sectoral_synergy) <- "cross-sectoral_synergy"
levels(cross_sectoral_synergy) <- data.frame(
  id = 0:15,
  trajectory = c("none", "climate", "economic", "economic + climate", "ecosystem",
                 "ecosystem + climate", "ecosystem + economic",
                 "ecosystem + economic + climate", "water", "water + climate",
                 "water + economic", "water + economic + climate",
                 "water + ecosystem", "water + ecosystem + climate",
                 "water + ecosystem + economic", "all sectors")
)

# plot ----
plot(c(intervention_priority, cross_sectoral_synergy))


# 11. Create bitfield registry ----
# socio_registry <- bf_registry(
#   name = "socioeconomic_planning",
#   description = "Socioeconomic analysis and intervention planning metrics"
# )
#
# # Add market distortion type (3 bits)
# socio_registry <- bf_test(
#   operator = "integer",
#   data = .rast(distortion_type),
#   x = distortion_type,
#   fields = list(sign = 0, exponent = 0, mantissa = 3),
#   registry = socio_registry
# )
#
# # Add market distortion magnitude (3 bits)
# socio_registry <- bf_test(
#   operator = "integer",
#   data = .rast(distortion_magnitude),
#   x = distortion_magnitude,
#   fields = list(sign = 0, exponent = 0, mantissa = 3),
#   registry = socio_registry
# )
#
# # Add GDP adjustment estimates (6 bits)
# socio_registry <- bf_test(
#   operator = "integer",
#   data = .rast(gdp_adjustment),
#   x = gdp_adjustment,
#   fields = list(sign = 0, exponent = 0, mantissa = 6),
#   registry = socio_registry
# )
#
# # Add Economic-Ecological Misalignment Index (8 bits)
# socio_registry <- bf_test(
#   operator = "integer",
#   data = .rast(misalignment),
#   x = economic_ecological_misalignment,
#   fields = list(sign = 0, exponent = 0, mantissa = 8),
#   registry = socio_registry
# )
#
# # Add system trajectory (3 bits)
# socio_registry <- bf_test(
#   operator = "integer",
#   data = .rast(system_trajectory),
#   x = system_trajectory,
#   fields = list(sign = 0, exponent = 0, mantissa = 3),
#   registry = socio_registry
# )
#
# # Add intervention priority (5 bits)
# socio_registry <- bf_test(
#   operator = "integer",
#   data = .rast(intervention_priority),
#   x = intervention_priority,
#   fields = list(sign = 0, exponent = 0, mantissa = 5),
#   registry = socio_registry
# )
#
# # Add cross-sectoral synergy (4 bits)
# socio_registry <- bf_test(
#   operator = "integer",
#   data = .rast(cross_sectoral_synergy),
#   x = cross_sectoral_synergy,
#   fields = list(sign = 0, exponent = 0, mantissa = 4),
#   description = c(
#     "Bitmap encoding of cross-sectoral synergies where:",
#     "- bit 3 (value 8): Water security synergy potential",
#     "- bit 2 (value 4): Ecosystem conservation synergy potential",
#     "- bit 1 (value 2): Economic development synergy potential",
#     "- bit 0 (value 1): Climate resilience synergy potential"
#   ),
#   registry = socio_registry
# )
#
#
# # Encode bitfield
# socio_field <- bf_encode(registry = socio_registry)
#
# # Create a raster version of the bitfield
# socio_rast <- rast(carrying_capacity)
# values(socio_rast) <- socio_field[,1]  # Assuming one column output
# names(socio_rast) <- "socioeconomic_planning_bitfield"
#
# # Plot the bitfield
# plot(socio_rast)

# 12. create plot items ----
#
dir.create(path = paste0(getwd(), "/figures/"))
toPlot <- c(mean_temp, annual_precip, soil_pH, soil_carbon, soil_clay,
            potential_npp, land_factor, soil_factor, slope_factor,
            carrying_capacity, animals_mean, carrying_capacity, exceedance_risk,
            exceedance, water_deficit, forage_deficit, soil_deficit,
            deficit_type, deficit_magnitude)

map(as.list(toPlot), .export_svg)
