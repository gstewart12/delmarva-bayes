
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLN", # three letter site code
  dir = "data", # where all data can be found
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

library(lubridate)
library(tidyverse)
# Packages required but not loaded: sf, smoothr, raster, stars

### Inputs

# Set working directory
wd <- file.path("data", settings$site, "processing")

paths <- list(
  #spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  point = "data/survey/flux_pts",
  # Path to site area
  delin = file.path(dirname(dirname(wd)), "site_info", "delineation"),
  # Path to site DEM
  dem = file.path(
    "data", "lidar", "DEM", "jackson_2007_1m", "2007_1m_DEM_filter.tif"
  ),
  # Path to half-hourly water level
  wtd = file.path(
    "~/Projects/Flux/inundation/data/output",
    paste0(str_replace(settings$site, "JL", "F"), "_waterlevel_dd.csv")
  ),
  # Path to point water level measurements
  wl = file.path(settings$dir, "Gas", "water_level.csv"),
  # Path to output files
  out = file.path(wd, "08_relative_elevation")
)
file_delin <- file.path(
  "~/Projects/Flux/spatial/data", settings$site,
  "processing/07_classification/classes_manual"
)

# Load metadata file
md <- yaml::read_yaml(file.path("data", settings$site, "metadata.yml"))
source("metadata/spatial_ref.R")


### Load input data

# Read AOI
#delin <- sf::read_sf(paths$delin)
# TODO: convert patch map to standalone delin somewhere else
delin <- file_delin %>%
  sf::read_sf() %>%
  filter(!class %in% c("UP", "WL")) %>%
  summarize() %>%
  smoothr::fill_holes(20)

# Read points of interest
point <- paths$point %>%
  sf::read_sf() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  filter(site == settings$site, type != "rain_gauge")
well <- filter(point, type == "well", is.na(zone))

# Read DEM
dem <- raster::raster(paths$dem) # DEM

# Crop DEM to site delineation
dem_f <- dem %>%
  raster::crop(raster::extent(sf::st_bbox(sf::st_buffer(delin, 10))))
  # Smooth to avoid influence of random error
  #raster::focal(raster::focalWeight(., 1, type = "Gauss"))

# Import water level data
data_wtd <- paths$wtd %>%
  read_csv(show_col_types = FALSE) %>%
  # Filter only complete water years
  filter(between(date, ymd("2017-10-01"), ymd("2020-09-30")))

# Find maximum water level measurement
# - using 99th percentile instead of true max
# - if it rains when the wetland is already "full", there is a brief spike in
#   water level during which overland flow is probably occurring
(wtd_max <- quantile(data_wtd$wtd_f, 0.99, na.rm = TRUE, names = FALSE))

# Find lowest point within site area
min_elev <- dem_f %>%
  raster::mask(delin) %>%
  raster::cellStats(min, na.rm = TRUE)
min_elev

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
  quantile(0.99, na.rm = TRUE, names = FALSE)
max_elev

# Determine relative elevation of well
# - equal to the difference between high water and measured water level
(well_rel_elev <- max_elev - wtd_max)

# Water level relative to lowest point in basin
wtd_corr <- data_wtd$wtd_f + well_rel_elev

# Water level as absolute elevation
wtd_abs <- wtd_corr + min_elev

# Check maximum flooded area
# - check any water level by changing 'wtd_max'
ggplot() +
  stars::geom_stars(
    data = rel_elev %>%
      raster::cut(c(0, well_rel_elev + wtd_max), include.lowest = TRUE) %>%
      stars::st_as_stars()
  ) +
  geom_sf(aes(color = zone), data = point) +
  geom_sf_text(aes(label = type), data = point, nudge_x = -7, nudge_y = 5) +
  scale_color_discrete(na.value = "black") +
  theme_void()

# Check flooded area at different water levels
ggplot() +
  stars::geom_stars(
    data = data_wtd %>%
      pull(wtd_f) %>%
      quantile() %>%
      map(
        ~ raster::cut(rel_elev, c(0, well_rel_elev + .x), include.lowest = TRUE)
      ) %>% 
      map(stars::st_as_stars) %>%
      lift_dl(c)(along = "wl")) +
  facet_wrap(~ wl) +
  theme_void()

# Write relative elevation raster to file
raster::writeRaster(
  rel_elev, file.path(paths$out, "rel_elev.tif"), overwrite = TRUE
)

# Write well elevation to file
# - good to also save it as a metadata entry


