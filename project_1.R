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

lc_meta <- data.frame(id = c(1, 2, 4:9),
                      var = c("trees", "grassland", "shrubs", "cropland", "built",
                              "bare", "water", "wetland"),
                      name = c("Forest", "Grassland", "Shrubland", "Cropland",
                               "Urban", "Bare", "Water", "Wetland"),
                      col = c("darkgreen", "yellowgreen", "wheat", "darkgoldenrod",
                              "firebrick2", "gray", "deepskyblue", "cyan2"))


# 1. get (bio)geophysical basis ----
#
## download and crop land cover data
landcover <- map(lc_meta$var, function(type){

  landcover(var = type, path = tempdir()) |>
    crop(area_extent) |> unlist()

}) |> rast()

## find dominant class
lc_dom <- which.max(landcover)
names(lc_dom) <- "landcover_dominating"

## get topography variables
elevation <- elevation_30s(country = "DEU", mask = FALSE, path = tempdir()) |>
  crop(area_extent)
names(elevation) <- "elevation"

slope <- terrain(elevation, "slope", unit = "degrees")

## get footprint
footprint <- footprint(year = 2009, path = tempdir()) |>
  crop(area_extent)
names(footprint) <- "human_footprint"

## distinguish extensive and intensive grassland (pasture)
lc_meta <- lc_meta |>
  rbind(data.frame(id = 3, var = "pasture", name = "Pasture", col = "green4"))
lc_meta <- lc_meta[order(lc_meta$id),]

lc_dom[lc_dom == 2 & footprint > 20] <- 9

lc_dom <- classify(x = lc_dom,
                   rcl = data.frame(is = c(1:9),
                                    becomes = c(1, 2, 4, 5, 6, 7, 8, 9, 3)))

levels(lc_dom) <- data.frame(id = lc_meta$id, landcover = lc_meta$name)
coltab(lc_dom) <- data.frame(value = lc_meta$id, col = lc_meta$col)

## plot
plot(c(elevation, slope, footprint, lc_dom))


# 2. generate simulated livestock densities ----
#
# this script generates simulated data matching the format of the GPW repository
# (https://zenodo.org/records/14933660).
#
## calculate factor affected by land
land_factor <-
  2.5 * landcover$grassland +
  1.0 * landcover$cropland +
  0.2 * landcover$trees
land_factor[landcover$built > 0.5] <- NA
land_factor[landcover$water > 0.5] <- NA

## calculate factor affected by elevation
# livestock production decreases with elevation, linear decline starting at 800m
elev <- 1.0 - pmax(0, values(elevation) - 800) / 1200
elev_factor <- elevation
values(elev_factor) <- elev

## calculate factor affected by slope
# steep slopes reduce grazing capacity, more than 20° slope significantly limits
# cattle production
slope_factor <- 1.0 - (slope / 45)
slope_factor[slope_factor < 0.1] <- 0.1  # Minimum factor on steep slopes

## calculate management factor
# random modification that is supposed to be due to management intensity
management_factor <- elevation
values(management_factor) <- runif(ncell(management_factor), 0.7, 1.2)

## calculate uncertainty
# uncertainty is assumed to be higher in mixed land use areas
land_mix <-
  (landcover$grassland * landcover$cropland) +
  (landcover$grassland * landcover$trees) +
  (landcover$cropland * landcover$trees)

# uncertainty model
uncertainty <- 0.15 +              # base uncertainty
  0.4 * land_mix +                 # land use heterogeneity
  0.01 * slope +                   # terrain complexity
  0.2 * management_factor          # management intensity

## calculate mean livestock density
animals_mean <- land_factor * elev_factor * slope_factor
names(animals_mean) <- "animals_mean"

## calculate standard deviation
animals_sd <- animals_mean * uncertainty
names(animals_sd) <- "animals_sd"

## calculate prediction interval bounds
# ±1.96 standard deviations
animals_lower <- animals_mean - 1.96 * animals_sd
names(animals_lower) <- "animals_lower"
animals_upper <- animals_mean + 1.96 * animals_sd
names(animals_upper) <- "animals_upper"

animals_lower[animals_lower < 0] <- 0

## plot
plot(c(landcover$grassland, animals_mean, animals_sd, animals_lower, animals_upper))


# 3. generate simulated ensemble runs from simulated layers ----
#
# this snippet takes the above simulated GPW mean and prediction interval and
# creates virtual but realistic ensemble model runs that would have generated
# these values.
model_types <- list(
  rf = list(name = "RandomForest", n_runs = 4, bias = 0.98, noise = 0.9),
  brt = list(name = "BoostedTrees", n_runs = 3, bias = 1.05, noise = 0.8),
  gam = list(name = "GAM", n_runs = 3, bias = 0.95, noise = 1.1),
  cnn = list(name = "CNN", n_runs = 2, bias = 1.02, noise = 0.7)
)

## simulate runs for each model type
ensemble <- map(names(model_types), function(name){

  thisModel <- model_types[[name]]

  temp <- map(1:thisModel$n_runs, function(ix){

    # create model-specific perturbation
    pertb <- animals_mean
    normNoise <- rnorm(ncell(animals_mean), mean = 0, sd = 1)

    # transform to match mean and standard deviation
    # For each pixel: mean*bias + (normNoise * noise_factor * animals_sd)
    perturbed_values <- thisModel$bias * values(animals_mean) + normNoise *
      thisModel$noise * values(animals_sd)
    values(pertb) <- perturbed_values
    pertb[pertb < 0] <- 0
    names(pertb) <- paste0(name, "_run", ix)
    return(pertb)

  }) |> rast()

}) |> rast()

# plot
plot(ensemble)


# 4. derive additional layers that are supposed to be reported ----
#
## calculate median
animals_median <- app(ensemble, fun = median, na.rm = TRUE)
names(animals_median) <- "animals_median"


## determine distribution type ----
dist_values <- app(ensemble, fun = function(x) {
  .get_distributionType(x)
})
names(dist_values) <- "dist_type"
levels(dist_values) <- data.frame(id = 0:7,
                                  distribution = c("normal", "lognormal", "beta",
                                                   "gamma", "weibull", "poisson",
                                                   "binomial", "other"))


## determine skewness and kurtosis ----
skewness <- app(ensemble, fun = function(x) {
  .get_skewness(x)
})
names(skewness) <- "skewness"
levels(skewness) <- data.frame(id = 0:7,
                               skewness = c("highly negative", "moder. negative",
                                            "slightly negative", "approx. symmetric",
                                            "slightly positive", "moder. positive",
                                            "highly positive", "other"))

kurtosis <- app(ensemble, fun = function(x) {
  .get_kurtosis(x)
})
names(kurtosis) <- "kurtosis"
levels(kurtosis) <- data.frame(id = 0:7,
                               kurtosis = c("flat", "moder. flat",
                                            "slightly flat", "mesokurtic",
                                            "slightly narrow", "moder. narrow",
                                            "narrow", "other"))

## simulate uncertainty source ----
# base on assumptions:
# 1. that areas with low footprint have worse data cover and thus more
#    uncertainty
# 2. that topographically more complex areas have more uncertainty
# 3. that pixels with more landcover classes have more uncertainty
# 4. that pixels with different landcover in the neighbourhood have more
#    uncertainty
uncertainty <- .simulate_uncertainty(landcover, elevation, footprint)

## determine model selection ----
model_names <- names(model_types)
model_avg <- list()
for (model in model_names) {
  layers <- grep(paste0("^", model), names(ensemble), value = TRUE)

  if (length(layers) > 0) {
    # calculate mean for this model type
    model_avg[[model]] <- app(subset(ensemble, layers), mean, na.rm = TRUE)
  }
}

model_result <- as.data.frame(c(rast(model_avg), animals_mean), xy = TRUE) |>
  rowwise() |>
  mutate(
    # calculate differences between each model and mean
    rf_diff = abs(rf - `animals_mean`),
    brt_diff = abs(brt - `animals_mean`),
    gam_diff = abs(gam - `animals_mean`),
    cnn_diff = abs(cnn - `animals_mean`),

    # determine which has minimum difference (0-based index)
    selection = which.min(c(rf_diff, brt_diff, gam_diff, cnn_diff)) - 1
  ) |>
  select(x, y, selection) |>
  as.matrix()

model_selection <- rasterize(x = model_result[,1:2], y = animals_mean, values = model_result[,3])
names(model_selection) <- "model_selection"
levels(model_selection) <- data.frame(id = 0:3,
                                      model = c("RF", "BRT", "GAM", "CNN"))

## determine model agreement ----
model_agreement <- app(rast(model_avg), function(x) {
  if (all(is.na(x))) return(NA)
  if (mean(x, na.rm = TRUE) == 0) return(0)
  1 - (sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE))  # 1 - CV = agreement
})

# transform to 8 categories (3 bits)
agreement_values <- values(model_agreement)
values(model_agreement) <- cut(agreement_values,
                               breaks = c(-Inf, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, Inf),
                               labels = FALSE) - 1
names(model_agreement) <- "model_agreement"
levels(model_agreement) <- data.frame(id = 0:7,
                                      agreement = c("severe disagr.", "high disagr.",
                                                    "moderate disagr.", "mild disagr.", "fair agr.",
                                                    "good agr.", "strong agr.", "perfect agr."))


## plot
plot(c(dist_values, model_selection, model_agreement, skewness, kurtosis, uncertainty$type, uncertainty$level))


# 5. Create bitfield registry ----
lstReg <- bf_registry(
  name = "livestock_density_modeling",
  description = "Bitfield encoding livestock density estimates with uncertainty metrics and model provenance.")

lstReg <- bf_map(
  protocol = "numeric",
  data = animals_median,
  registry = lstReg,
  name = "Median livestock density",
  x = animals_median,
  fields = list(exponent = 2, significand = 5),
  na.val = 0)

lstReg <- bf_map(
  protocol = "integer",
  data = animals_sd,
  registry = lstReg,
  name = "Standard deviation",
  x = animals_sd,
  fields = list(significand = 7),
  na.val = 0)

lstReg <- bf_map(
  protocol = "category",
  data = dist_values,
  registry = lstReg,
  name = "Distribution type",
  x = distribution,
  na.val = 7)

lstReg <- bf_map(
  protocol = "category",
  data = skewness,
  registry = lstReg,
  name = "Skewness",
  x = skewness,
  na.val = 7)

lstReg <- bf_map(
  protocol = "category",
  data = kurtosis,
  registry = lstReg,
  name = "Kurtosis",
  x = kurtosis,
  na.val = 7)

lstReg <- bf_map(
  protocol = "integer",
  data = uncertainty,
  registry = lstReg,
  name = "Confidence level",
  x = level,
  fields = list(significand = 4))

lstReg <- bf_map(
  protocol = "category",
  data = uncertainty,
  registry = lstReg,
  name = "Uncertainty source",
  x = type)

lstReg <- bf_map(
  protocol = "category",
  data = model_selection,
  registry = lstReg,
  name = "Model selection",
  x = model,
  na.val = 0)

lstReg <- bf_map(
  protocol = "category",
  data = model_agreement,
  registry = lstReg,
  name = "Model Agreement Index",
  x = agreement,
  na.val = 0)

field <- bf_encode(registry = lstReg)
rst_field <- rast(lc_dom, vals = field, nlyrs = 2, names = names(field))


# 5. export items ----
#
dir.create(path = paste0(getwd(), "/figures/"))
toPlot <- c(elevation, slope, footprint, lc_dom, landcover$grassland,
            animals_mean, animals_median, animals_sd, animals_lower, animals_upper,
            ensemble$rf_run1, ensemble$brt_run1, ensemble$gam_run1,
            ensemble$cnn_run1, dist_values, model_selection, model_agreement,
            skewness, kurtosis, uncertainty$type, uncertainty$level)

map(as.list(toPlot), .export_svg)

writeRaster(x = rst_field, filename = "lstBitfield.tif", overwrite = TRUE)
writeRaster(x = animals_mean, filename = "animals_mean.tif", overwrite = TRUE)

bf_export(registry = lstReg, format = "yaml", file = paste0(getwd(), "/meta/lstMeta.yml"))
saveRDS(object = lstReg, file = paste0(getwd(), "/lstRegistry.rds"))
