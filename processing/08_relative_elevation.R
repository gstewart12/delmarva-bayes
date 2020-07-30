
rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  dir = "/Users/Graham/Desktop/DATA", # where all data can be found
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

library(tidyverse)


### Inputs

# Set working directory
wd <- file.path(
  "/Users/Graham/Desktop/DATA/Flux", settings$site, "analysis", "spatial"
)

paths <- list(
  spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  point = "/Users/Graham/Desktop/DATA/Spatial/survey/flux_pts",
  # Path to site area
  delin = file.path(dirname(dirname(wd)), "site_info", "delineation"),
  # Path to site DEM
  # - can also use DEM from 01_clip_imagery but this one is already on site grid
  dem = file.path(settings$dir, "Spatial", "DEM", "1m_Relief_invert.tif"),
  # Path to reference grid
  grid = file.path(wd, "03_feature_extraction", "lidar", "dem.tif"),
  # Path to half-hourly water level
  wtd = file.path(
    "/Users/Graham/Desktop/DATA/PT/output", settings$site, 
    paste0(settings$site, "_all.csv")
  ),
  # Path to point water level measurements
  wl = file.path(settings$dir, "Gas", "water_level.csv"),
  # Path to output files
  out = file.path(wd, "08_relative_elevation")
)

# Load reference files
source(paths$spatial_ref)

raster_rescale <- function(x, ...) {
  values <- raster::values(x)
  new_values <- scales::rescale(values, ...)
  raster::setValues(x, new_values)
}


### Load input data

# Read AOI
delin <- sf::read_sf(paths$delin)

# Read points of interest
point <- paths$point %>%
  sf::read_sf() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::filter(site == settings$site, type != "rain_gauge")
well <- dplyr::filter(point, type == "well", is.na(zone))

# Read DEM
dem <- raster::raster(paths$dem) # DEM

# Smooth and resample DEM to site grid
grid <- raster::raster(paths$grid)
dem_f <- dem %>%
  raster::crop(raster::extent(sf::st_buffer(delin, 10))) %>%
  # Smooth to avoid influence of random error
  #raster::focal(raster::focalWeight(., 1, type = "Gauss")) %>%
  raster::resample(grid, method = "bilinear")

# Import water level data
data_wtd <- paths$wtd %>%
  readr::read_csv(
    guess_max = 6000, col_types = readr::cols(.default = readr::col_guess())
  ) %>%
  #dplyr::select(timestamp, wtd_f, dry_flag, p) %>%
  # Sub-daily measurements are not necessary here
  dplyr::mutate(date = lubridate::date(timestamp)) %>% 
  dplyr::group_by(date) %>% 
  dplyr::summarize(wtd = mean(wtd_f, na.rm = TRUE), .groups = "drop")


# Find median water level measurement
# - this just needs to be somewhere between completely dry and wetlands merging
#wtd_med <- median(data_wtd$wtd_f, na.rm = TRUE)

# Find maximum water level measurement
# - using 99th percentile instead of true max
# - if it rains when the wetland is already "full", there is a brief spike in
#   water level during which overland flow is probably occurring
wtd_max <- quantile(data_wtd$wtd, 0.99, na.rm = TRUE, names = FALSE)
wtd_max

# Find lowest point within site area
min_elev <- dem_f %>%
  raster::mask(delin) %>%
  raster::cellStats(min, na.rm = TRUE)
min_elev
#well_elev <- raster::extract(dem, well)

# min_elev <- dem_aoi %>%
#   raster::cut(c(0, well_elev + wtd_med), include.lowest = TRUE) %>%
#   #raster::reclassify(c(-Inf, well_elev + wtd_max, 1, well_elev + wtd_max, Inf, 0)) %>% 
#   stars::st_as_stars() %>%
#   sf::st_as_sf(as_points = FALSE, merge = TRUE) %>%
#   tibble::as_tibble() %>% 
#   dplyr::filter(sf::st_intersects(geometry, well, sparse = FALSE)) %>%
#   sf::st_as_sf() %>%
#   raster::mask(dem_aoi, .) %>%
#   raster::cellStats(min, na.rm = TRUE)

# Normalize DEM to minimum point
rel_elev <- dem_f - min_elev
raster::plot(rel_elev)

# Find highest relative elevation within site area
# - high water level is when 99% of area is flooded
# - this allows for some wetland area that is only rarely flooded 
#   (but still has wetland soils, vegetation, etc.) 
max_elev <- rel_elev %>%
  raster::mask(delin) %>%
  raster::values() %>%
  #max(na.rm = TRUE)
  quantile(0.99, na.rm = TRUE, names = FALSE)
max_elev

# Determine relative elevation of well
# - equal to the difference between high water and measured water level
(well_rel_elev <- max_elev - wtd_max)

# Check maximum flooded area
# - check flooding extent of any water level by changing 'wtd_max'
ggplot2::ggplot() +
  stars::geom_stars(
    data = rel_elev %>%
      raster::cut(c(0, well_rel_elev + wtd_max), include.lowest = TRUE) %>%
      stars::st_as_stars()
  ) +
  ggplot2::geom_sf(ggplot2::aes(color = zone), data = point) +
  ggplot2::geom_sf_text(
    ggplot2::aes(label = type), data = point, nudge_x = -7, nudge_y = 5
  ) +
  ggplot2::theme_void()


# Write relative elevation raster to file
raster::writeRaster(
  rel_elev, file.path(paths$out, "rel_elev.tif"), overwrite = TRUE
)

# Write well elevation to file
# - good to also save it as a metadata entry



# Read water level point measurement data
data_wl <- readr::read_csv(paths$wl) %>%
  dplyr::mutate(site = stringr::str_replace_all(site, "F", "JL")) %>%
  dplyr::filter(site == settings$site) %>%
  dplyr::select(date, site, zone, wl_ch, wtd_zone, flooded) %>%
  # Add site-level WTD
  dplyr::left_join(data_wtd, by = "date")
data_wl %>% print(n = 50)



