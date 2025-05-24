library(terra)
library(geodata)
library(tidyverse)
library(bitfield)

source("functions.R")

set.seed(42)

xmin <- 7.6
xmax <- 9.3
ymin <- 47.5
ymax <- 48.8
area_extent <- ext(xmin, xmax, ymin, ymax)

# 1. Get (bio)geophysical basis ----
#
## Download and crop climate, soil and elevation data ----
tavg <- worldclim_global(var = "tavg", res = 0.5, path = tempdir()) |>
  crop(area_extent)
tmax <- worldclim_global(var = "tmax", res = 0.5, path = tempdir()) |>
  crop(area_extent)
prec <- worldclim_global(var = "prec", res = 0.5, path = tempdir()) |>
  crop(area_extent)
soil_pH <- soil_world(var = "phh2o", depth = 5, path = tempdir()) |>
  crop(area_extent)
soil_carbon <- soil_world(var = "soc", depth = 5, path = tempdir()) |>
  crop(area_extent)
soil_clay <- soil_world(var = "clay", depth = 15, path = tempdir()) |> # for some reason the depth = 5 layer doesn't work currently
  crop(area_extent)

elevation <- elevation_30s(country = "DEU", mask = FALSE, path = tempdir()) |>
  crop(area_extent)
names(elevation) <- "elevation"
slope <- terrain(elevation, "slope", unit = "degrees")
aspect <- terrain(elevation, "aspect", unit = "degrees")

# Summarize climate variables
mean_temp <- mean(tavg)
names(mean_temp) <- "yearly_average_temperature"
max_month <- max(tmax)
names(max_month) <- "warmest_month"
annual_precip <- sum(prec)
names(annual_precip) <- "yearly_total_precipitation"

# Plot
plot(c(mean_temp, annual_precip, soil_pH, soil_carbon, soil_clay, ))


# 2. Decode bitfield from project_1 to retrieve layers ----
# Assuming these layers were decoded from Project 1's bitfield
# and are available in the environment
# animals_mean
# animals_sd
# dist_values
# skewness
# kurtosis
# landcover


# 3. calculate carrying capacity ----
#
# calculate NPP using Miami model with climate data
npp_temp <- 3000 / (1 + exp(1.315 - 0.119 * mean_temp))
npp_precip <- 3000 * (1 - exp(-0.000664 * annual_precip))
potential_npp <- min(npp_temp, npp_precip)
names(potential_npp) <- "potential_npp"

# determine landcover factor from simplified classifications
land_factor <-
  1.5 * landcover$grassland +                       # grassland most productive for livestock
  1.0 * landcover$cropland +                        # cropland provides some grazing
  0.7 * landcover$trees +                           # forest has limited grazing
  0.5 * landcover$shrubs +                          # shrubland has minimal forage
  0.2 * landcover$bare                              # bare land can support minimal vegetation
names(land_factor) <- "land_support"

# determine soil factor
pH_factor <- 1 - (0.2 * abs(soil_pH - 6.5) / 2)     # optimal pH for grassland is around 6-7, with decreasing productivity away from this range
carbon_factor <- 0.7 + (0.3 * (soil_carbon / 200))  # organic carbon improves soil fertility and moisture retention
clay_factor <- 1 - (0.2 * abs(soil_clay - 25) / 20) # optimal clay is around 20-30%, with decreased productivity at extremes
soil_factor <- (pH_factor + carbon_factor + clay_factor) / 3
names(soil_factor) <- "soil_support"

# determine slope factor (Otieno, 2019)
slope_factor <- 1.0 - (min(slope, 45) / 45)^1.5
names(slope_factor) <- "slope_limitation"

# define sustainable utilization rate
proper_use_factor <- 0.25

# Calculate actual carrying capacity
usable_npp <- potential_npp *
  land_factor *
  soil_factor *
  slope_factor *
  proper_use_factor

# Convert NPP (g/m²/yr) to livestock units per hectare assuming 2300 kg dry matter per livestock unit per year
carrying_capacity <- usable_npp * 10 / 2300
names(carrying_capacity) <- "carrying_capacity"

# plot
plot(c(potential_npp, land_factor, soil_factor, slope_factor, carrying_capacity))


# 4. Determine carrying capacity exceedance ----
exceedance_prob <- .exceedance_probability(animals_mean, animals_sd, dist_values,
                                           skewness, kurtosis, carrying_capacity)

exceedance_risk <- app(exceedance_prob, function(x) {
  round(x * 63)                                     # scale 0-1 to 0-63
})
names(exceedance_risk) <- "exceedance_risk"

exceedance <- animals_mean / carrying_capacity
names(exceedance) <- "exceedance_magnitude"

# Plot
plot(c(animals_mean, carrying_capacity, exceedance_risk, exceedance))


# 5. Calculate resource limitations ----
## Water limitation ----
water_demand <- animals_mean * 45                   # water need per livestock unit (liters/day)

# south-facing slopes have increased evaporation
south_factor <- (1 - abs(aspect - 180) / 180) * 0.3 + 0.7  # 0.7-1.0 range (lower = south-facing)

# get the driest month from the precipitation data
driest_month <- min(prec) * 0.7                    # find minimum across all 12 months
driest_daily <- driest_month / 30                    # convert to daily precipitation

# increased animal water requirements during heat periods
heat_stress_factor <- 1 + ((max_month - 15) / 4)

# apply to demand
water_demand <- water_demand * heat_stress_factor

# calculate water availability - baseline calculation
total_water <- elevation                             # template raster
values(total_water) <- values(driest_daily) * 0.4 * 10000  # precipitation with retention factor

# calculate animal-available water (only a small fraction of total)
available_water <- total_water * 0.05         # ~5% becomes available for drinking

# areas with cropland have reduced water availability due to competition
crop_intensity <- (landcover$cropland)^1.5           # non-linear relationship with intensity
crop_extraction <- 1 - (crop_intensity * 0.3)        # up to 30% reduction in high-intensity cropland

# apply all modifying factors to animal-available water
adjusted_water <- available_water *
  (1 - 0.8 * (values(slope) / 45)) *                 # topographic slope effect
  values(south_factor) *                             # south-facing slope effect
  values(crop_extraction)                            # agricultural water extraction

# calculate deficit ratio using heat-adjusted demand
water_deficit <- water_demand / adjusted_water
water_deficit[water_deficit < 1] <- 0                # no deficit where supply exceeds demand
names(water_deficit) <- "water_deficit"


## Forage limitation ----
forage_demand <- animals_mean * 10                  # kg dry matter per livestock unit per day

forage_availability <- potential_npp * 0.0005 *     # kg/m²/yr
  0.7 * landcover$grassland +
  0.5 * landcover$cropland +
  0.2 * landcover$trees
forage_availability <- forage_availability * 10000  # annual production per ha
forage_deficit <- (forage_demand * 365) / (forage_availability + 0.001)  # add small constant to avoid division by zero
forage_deficit[forage_deficit < 1] <- 0
names(forage_deficit) <- "forage_deficit"

## Soil fertility limitation ----
soil_ph_factor <- 1 - abs(soil_pH - 6.5) / 3.5      # optimal pH around 6.5
soil_carbon_factor <- soil_carbon / 200             # higher carbon is better for soil resilience
soil_clay_factor <- 1 - abs(soil_clay - 25) / 75    # optimal clay around 25%

soil_fertility <- min(c(soil_ph_factor, soil_carbon_factor, soil_clay_factor))
soil_fertility <- ifel(soil_fertility < 0, 0, soil_fertility)  # ensure non-negative

base_impact <- 0.2
slope_multiplier <- 1 + (slope / 45)^1.5            # steeper slope means more erosion risk
soil_demand <- animals_mean * base_impact * slope_multiplier

soil_deficit <- soil_demand / (soil_fertility + 0.001)  # add small constant to avoid division by zero
soil_deficit[soil_deficit < 1] <- 0
names(soil_deficit) <- "soil_deficit"

## Resource limitation type ----
# 0=no limitation, 1=water, 2=forage, 3=soil
deficit_type <- ifel(water_deficit > forage_deficit & water_deficit > soil_deficit & water_deficit > 0, 1,
                      ifel(forage_deficit > soil_deficit & forage_deficit > 0, 2,
                           ifel(soil_deficit > 0, 3, 0)))
names(deficit_type) <- "deficit_type"
levels(deficit_type) <- data.frame(id = 0:3,
                                    limitation = c("none", "water", "forage", "soil"))

## Resource limitation magnitude ----
deficit_magnitude <- ifel(deficit_type == 1, log10(water_deficit + 1),
                           ifel(deficit_type == 2, log10(forage_deficit + 1),
                                ifel(deficit_type == 3, log10(soil_deficit + 1), 0)))

# Scale to 0-7 range for 3-bit encoding
deficit_magnitude <- (deficit_magnitude / max(values(deficit_magnitude), na.rm = TRUE)) * 7
names(deficit_magnitude) <- "deficit_magnitude"

# Plot
plot(c(water_deficit, forage_deficit, soil_deficit, deficit_type, deficit_magnitude))


# 6. Create bitfield ----
#
this also needs an update of the previous bitflag for distribution types.
# # Create a registry for the ecological analysis
# eco_registry <- bf_registry(
#   name = "ecological_analysis",
#   description = "Ecological carrying capacity analysis derived from livestock density data"
# )
#
# # Add carrying capacity (7 bits)
# eco_registry <- bf_test(
#   operator = "numeric",
#   data = .rast(carrying_capacity),
#   x = carrying_capacity,
#   fields = list(sign = 0, exponent = 0, mantissa = 7),
#   registry = eco_registry
# )
#
# # Add carrying capacity exceedance risk (6 bits)
# eco_registry <- bf_test(
#   operator = "numeric",
#   data = .rast(exceedance_risk),
#   x = exceedance_risk,
#   fields = list(sign = 0, exponent = 0, mantissa = 6),
#   registry = eco_registry
# )
#
# # Add carrying capacity exceedance magnitude (6 bits)
# eco_registry <- bf_test(
#   operator = "numeric",
#   data = .rast(exceedance),
#   x = exceedance,
#   fields = list(sign = 0, exponent = 0, mantissa = 6),
#   registry = eco_registry
# )
#
# # Add resource limitation type (2 bits)
# eco_registry <- bf_test(
#   operator = "integer",
#   data = .rast(deficit_type),
#   x = resource_limitation_type,
#   fields = list(sign = 0, exponent = 0, mantissa = 2),
#   registry = eco_registry
# )
#
# # Add resource limitation magnitude (3 bits)
# eco_registry <- bf_test(
#   operator = "integer",
#   data = .rast(deficit_magnitude),
#   x = resource_limitation_magnitude,
#   fields = list(sign = 0, exponent = 0, mantissa = 3),
#   registry = eco_registry
# )
#
# # Encode the bitfield
# eco_field <- bf_encode(registry = eco_registry)
#
# # Create a raster version of the bitfield
# eco_rast <- rast(carrying_capacity)
# values(eco_rast) <- eco_field[,1]  # Assuming one column output
# names(eco_rast) <- "ecological_analysis_bitfield"
#
# # Plot the bitfield
# plot(eco_rast)

# 7. create plot items ----
#
dir.create(path = paste0(getwd(), "/figures/"))
toPlot <- c(mean_temp, annual_precip, soil_pH, soil_carbon, soil_clay,
            potential_npp, land_factor, soil_factor, slope_factor,
            carrying_capacity, animals_mean, carrying_capacity, exceedance_risk,
            exceedance, water_deficit, forage_deficit, soil_deficit,
            deficit_type, deficit_magnitude)

map(as.list(toPlot), .export_svg)
