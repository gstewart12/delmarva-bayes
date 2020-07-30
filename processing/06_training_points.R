

settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLN", # three letter site code
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

# References
# 
# Millard, K., & Richardson, M. (2015). On the Importance of Training Data 
# Sample Selection in Random Forest Image Classification: A Case Study in 
# Peatland Ecosystem Mapping. Remote Sensing, 7(7), 8489–8515. 
# https://doi.org/10.3390/rs70708489
# 
# Smith, A. (2010). Image segmentation scale parameter optimization and land 
# cover classification using the Random Forest algorithm. Journal of Spatial 
# Science, 55(1), 69–79. https://doi.org/10.1080/14498596.2010.487851

control <- list(
  sample_n_init = 1000,
  # Choosing sample size a priori rather than per area - just need to make sure
  #   there will be enough samples in each class for RF classifier
  sample_n_final = 200,
  # Minimum distance allowed between points
  buffer_dist = 5
)

library(sf)
library(tidyverse)

raster_rescale <- function(x, ...) {
  values <- raster::values(x)
  new_values <- scales::rescale(values, ...)
  raster::setValues(x, new_values)
}

geom_stars_rgb <- function(data, r = 1, g = 2, b = 3, ...) {
  data <- data %>%
    dplyr::select(
      r = dplyr::all_of(r), g = dplyr::all_of(g), b = dplyr::all_of(b)
    ) %>%
    dplyr::mutate(rgb = grDevices::rgb(r, g, b, maxColorValue = 255))
  stars::geom_stars(data = data, ggplot2::aes(x = x, y = y, fill = rgb), ...)
}


# Set working directory
wd <- file.path(
  "/Users/Graham/Desktop/DATA/Flux", settings$site, "analysis", "spatial"
)

paths <- list(
  spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  seeds = "/Users/Graham/Desktop/DATA/Flux/tools/reference/seeds.R",
  # Path to site area
  delin = file.path(wd, "02_site_delineation", "delineation"),
  point = "/Users/Graham/Desktop/DATA/Spatial/survey/flux_pts",
  # Output path
  out = file.path(wd, "06_training_points")
)

# Load reference files
source(paths$spatial_ref)
source(paths$seeds)
seeds <- purrr::pluck(seeds, settings$site, "training_points")


# File paths for all outputs
files_out <- list(
  sample = file.path(paths$out, "sample_points", "sample_points.shp"),
  ref = file.path(paths$out, "reference_points", "ref_points.shp"),
  train = file.path(paths$out, "training_points", "train_points.shp")
)


### Load input data


# Read AOI
delin <- sf::read_sf(paths$delin)

# Read points of interest
point <- sf::read_sf(paths$point)
tower <- point %>% 
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::filter(site == settings$site, type == "tower") %>%
  sf::st_geometry()


# Selection of training areas via random sampling (Millard & Richardson 2015)

# Perform initial random sample of points
set.seed(seeds$sample) # SET THE SEED
sample_init <- delin %>%
  # Don't sample too close to the edge - segments might extend outside here
  sf::st_buffer(-5) %>%
  sf::st_sample(control$sample_n_init) %>%
  sf::st_as_sf() %>%
  tibble::rowid_to_column(var = "ID")
  # Indicate which segment each point lies in
  #sf::st_join(segments_poly) %>%
  #dplyr::rename(segment = segments_init)

# Check sample
ggplot2::ggplot() +
  ggplot2::geom_sf(data = delin, col = "gray40", fill = NA, lwd = 1) +
  ggplot2::geom_sf(data = sample_init, pch = 3, col = "red", alpha = 0.7) +
  ggplot2::theme_void()

# Write points to file
sf::write_sf(sample_init, files_out$sample)

# Prune random points to a minimum separated distance
# - consider removing this requirement if using 1-per-segment rule anyway
# - although it still ensures that samples are somewhat spread out
# - keeping for now
sample_prune <- sample_init
i <- 1 # iterator start
repeat ({
  # Create buffer around i-th point, determine intersection
  reject <- sample_prune %>%
    sf::st_intersects(
      sf::st_buffer(dplyr::slice(sample_prune, i), control$buffer_dist)
    ) %>%
    tibble::as_tibble() %>%
    # Remove the first element (it is the original point)
    dplyr::slice(-1) %>% 
    # Take only the row.id column as a vector
    dplyr::pull(row.id) 
  
  # Exclude any points not observing social distancing
  if (!rlang::is_empty(reject)) {
    sample_prune <- dplyr::slice(sample_prune, -reject)
  } 
  
  if (i >= nrow(sample_prune)) {
    break # the end was reached; no more points to process
  } else {
    i <- i + 1 # repeat
  }
})

# Make sure there is only one point per segment
# - NO! better to keep training points independent from segmentation
# - not a problem if this results in slightly fewer training segments
#sample_prune <- sample_prune %>% 
#  dplyr::group_by(segment) %>% 
#  dplyr::slice(1) %>%
#  dplyr::ungroup()

# Check pruning
ggplot2::ggplot() +
  ggplot2::geom_sf(data = delin, col = "gray40", fill = NA, lwd = 1) +
  ggplot2::geom_sf(data = sample_prune, pch = 3, col = "red", alpha = 0.7) +
  ggplot2::theme_void()

# Subset to final sample size
set.seed(seeds$subset) # SET THE SEED
sample_points <- sample_prune %>% 
  tibble::as_tibble() %>% 
  # Optimize representation of samples near flux tower
  # - weight = 1 / distance^2
  dplyr::mutate(dist = 1 / as.numeric(sf::st_distance(x, tower))^2) %>% 
  # Final subset selection is weighted using inverse distance from tower
  dplyr::slice_sample(weight_by = dist, n = control$sample_n_final) %>%
  dplyr::select(-dist) %>%
  sf::st_as_sf()

# Check final sample
ggplot2::ggplot() +
  ggplot2::geom_sf(data = delin, col = "gray40", fill = NA, lwd = 1) +
  ggplot2::geom_sf(data = sample_points, pch = 3, col = "red", alpha = 0.7) +
  ggplot2::theme_void()

# Write final sample - one set for reference points, one for classification
sf::write_sf(sample_points, files_out$ref)
# Add class & certainty fields to training points
# - CAREFUL - don't overwrite existing file with manual classifications!
sf::write_sf(
  dplyr::mutate(
    sample_points, 
    class = NA_character_, 
    subclass = NA_character_,
    certainty = NA_integer_, 
    .after = ID
  ), 
  files_out$train
)
# Best to classify points in case segmentation method changes
# - segments can be used as reference in discerning the patch a point resides in


# Probably end script here


# Intersect points with segments
sample_segments <- segments_poly %>%
  sf::st_join(sample_points) %>%
  dplyr::select(-segments_init) %>%
  # Only keep segments that intersected with a sample point
  dplyr::filter(!is.na(ID))

# Write a set of bare training segments for reference
sf::st_write(sample_segments, files_out$ref_segm, append = FALSE)

# Check sample segments
ggplot2::ggplot() +
  geom_stars_rgb(data = img_bands, r = 3, g = 4, b = 5) +
  ggplot2::geom_sf(data = delin, col = "gray80", fill = NA, lwd = 1) +
  ggplot2::geom_sf(data = sample_segments, col = "gray80", fill = NA) +
  ggplot2::scale_fill_identity() +
  ggplot2::theme_void()

# Summarize data by segment to help with classification 
segment_bands <- img_bands %>% 
  # Treat cells as points so adjacent cells are not included
  sf::st_join(sample_segments, what = "inner", as_points = TRUE) %>%
  tibble::as_tibble() %>% 
  dplyr::select(-x, -y) %>% 
  dplyr::group_by(ID, segment, geometry) %>% 
  dplyr::summarize(
    dplyr::across(.fns = median), n = dplyr::n(), .groups = "drop"
  ) %>%
  # Rounding values makes a prettier attribute table in QGIS
  dplyr::mutate(dplyr::across(c(red:nir), ~ as.integer(round(.x, 0)))) %>%
  dplyr::mutate(dplyr::across(c(chm:dem), ~ round(.x, 2))) %>%
  sf::st_as_sf() %>%
  # Remove double-counted segments (if any)
  dplyr::distinct(segment, .keep_all = TRUE)

# Calculate average Euclidean distances for cells in each segment
segment_dists <- img_bands %>% 
  # Treat cells as points so adjacent cells are not included
  sf::st_join(sample_segments, what = "inner", as_points = TRUE) %>%
  tibble::as_tibble() %>% 
  dplyr::select(-segment, -x, -y, -geometry) %>% 
  # Rescale to balance distances across bands
  dplyr::mutate(dplyr::across(-ID, scales::rescale)) %>%
  tidyr::pivot_longer(-ID) %>% 
  dplyr::group_by(ID, name) %>%
  dplyr::mutate(cell = seq_along(value)) %>%
  dplyr::group_by(cell, .add = TRUE) %>% 
  dplyr::mutate(mean = mean(value), dist = abs(value - mean)^2) %>% 
  dplyr::group_by(ID, cell) %>% 
  dplyr::summarize(d = sqrt(sum(value)), .groups = "drop_last") %>% 
  dplyr::summarize(d = round(mean(d), 2), .groups = "drop")

# Join Euclidean distance with main data frame, prep for classification
training_segments <- segment_bands %>%
  dplyr::left_join(segment_dists, by = "ID") %>%
  # Add class & certainty fields
  dplyr::mutate(class = NA_character_, certainty = NA_integer_, .after = ID) %>%
  dplyr::mutate(dplyr::across(where(is.double), as.character))

# Write training segments for classification
sf::st_write(training_segments, files_out$train_segm)

# Check some summary statistics of training sample
sample_segments %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    area = sf::st_area(geometry),
    tower_dist = sf::st_distance(geometry, tower),
    edge_dist = sf::st_distance(geometry, sf::st_cast(delin, "MULTILINESTRING"))
  ) %>%
  dplyr::summarize(
    n = n(), 
    n_unique = dplyr::n_distinct(segment),
    p_tot_segm = n_unique / n_segm,
    avg_tower_dist = mean(tower_dist),
    avg_edge_dist = mean(edge_dist),
    avg_area = mean(area),
    tot_area = sum(area),
    p_site_area = tot_area / sf::st_area(delin)
  )

# Check sample segments
ggplot2::ggplot() +
  ggplot2::geom_sf(data = delin, col = "gray80", fill = NA, lwd = 1) +
  ggplot2::geom_sf(
    data = dplyr::left_join(segment_bands, segment_dists), aes(fill = d), 
    col = "gray80"
  ) +
  ggplot2::scale_fill_distiller(palette = "Spectral") +
  ggplot2::theme_void()

# Get VERY ROUGH initial estimation of class balance in training segments
segment_bands %>%
  tibble::as_tibble() %>%
  dplyr::mutate(
    class_est = dplyr::case_when(
      chm > 5 ~ "frs", 
      red > 80 ~ "veg",
      nir < 100 ~ "wtr",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(class_est)) %>%
  dplyr::group_by(class_est) %>%
  dplyr::summarize(p = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(p = p / sum(p))


# Finished - training sample classification is done manually in QGIS
