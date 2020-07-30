
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLA", # three letter site code
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

library(tidyverse)

# This script is a product of Nate Jones' work in:
# https://github.com/FloodHydrology/DEM_Inundate
# https://github.com/FloodHydrology/DMV_Spatial_Analysis

# Outline
# 1. Load DEM
# 2. Crop DEM to site AOI
# 3. Fill single cell pits in DEM
# 4. Smooth DEM - gaussian filter; sigma = 2 (or 3?)
# 5. Run stochastic depression analysis; rmse = 0.18, range = 10, iter = 500
# 6. Reclassify SDA as depression or upland; thresh. prob. = 0.1 (or 0.05?)

# Define relevant working directories
workspace <- "/Users/Graham/Desktop/DATA/Flux/scratch"
#wbt_path <- "/Users/Graham/WBT/whitebox_tools"
path_out <- file.path("/Users/Graham/Desktop/DATA/Flux", settings$site, "site_info")

# Load reference files
source("~/Desktop/DATA/Flux/tools/reference/seeds.R")
source("~/Desktop/DATA/Flux/tools/reference/site_metadata.R")
source("~/Desktop/DATA/Flux/tools/reference/spatial.R")

seeds <- purrr::pluck(seeds, settings$site)

# Download relevant data 
# DEMs: 2007_1m_DEM, 1m_Relief_invert
dem <- raster::raster("~/Desktop/DATA/Spatial/DEM/1m_Relief_invert.tif")
aoi <- sf::st_read(
  file.path("~/Desktop/DATA/Flux", settings$site, "site_info", "aoi")
)
aoi <- sf::st_transform(aoi, crs = sp::proj4string(dem))

dem_crop <- raster::crop(dem, aoi)
raster::plot(dem_crop)
raster::writeRaster(dem_crop, file.path(workspace, "dem.tif"), overwrite = TRUE)

# Prep DEM for Analysis ========================================================

# Fill single cell pits
whitebox::wbt_fill_single_cell_pits(
  dem = file.path(workspace, "dem.tif"), 
  output = file.path(workspace, "dem_fill.tif")
)

# Read raster back into workspace
dem_fill <- raster::raster(file.path(workspace, "dem_fill.tif"))

# Apply simple gausian filter to smooth random errors from DEM
dem_filter <- raster::focal(dem_fill, w = raster::focalWeight(dem, 3, "Gauss"))

# Export filtered DEM to workspace
raster::writeRaster(
  dem_filter, file.path(workspace, "dem_filter.tif"), overwrite = TRUE
)


# Identify depression in landscape =============================================

# Identify depressions using WhiteBox GAT Monte Carlo Approach

sda_params <- list(
  rmse = 0.18, # RMSE of 18 cm
  range = 10, 
  iter = 10000 # high number to ensure stable results
)

# Takes a couple minutes to run
set.seed(seeds$delineation$sda) # doesn't look like this ensures stable results
whitebox::wbt_stochastic_depression_analysis(
  dem = file.path(workspace, "dem_filter.tif"),
  output = file.path(workspace, "dem_sda.tif"),
  rmse = sda_params$rmse,
  range = sda_params$range,
  iterations = sda_params$iter
)

# Define depression based on threshold of occurence

# Reclassify raster 
dem_sda <- raster::raster(file.path(workspace, "dem_sda.tif"))
raster::plot(dem_sda)

t <- 0.1 # set cutoff threshold - more aimed at removing definite upland
rcl_matrix <- matrix(
  c(0, t, NA, t, 1, 1), 
  ncol = 3, byrow = TRUE
)
dem_reclass <- raster::reclassify(dem_sda, rcl_matrix, include.lowest = TRUE)
raster::plot(dem_reclass)


# Delineate wetland area =======================================================

# Easier to work with vectors instead of rasters at this point
polys <- dem_reclass %>%
  stars::st_as_stars() %>%
  sf::st_as_sf(as_points = FALSE, merge = TRUE)
plot(polys)

# The largest polygon is the basin
poly <- polys %>%
  dplyr::mutate(area = sf::st_area(geometry)) %>%
  dplyr::slice_max(area) %>%
  dplyr::select(-area)

# Clean up polygon
poly_fill <- smoothr::fill_holes(poly, 200)
plot(poly_fill)

# "Pinch off" peninsulas and isthmi with negative buffer
poly_trim <- poly_fill %>% 
  sf::st_buffer(-10, joinStyle = "MITRE") %>% 
  sf::st_cast("POLYGON") %>% 
  dplyr::mutate(area = sf::st_area(geometry)) %>% 
  dplyr::slice_max(area) %>% 
  dplyr::select(-area) %>% 
  # Apply positive buffer to recover lost area
  sf::st_buffer(10, joinStyle = "MITRE")
# This also helps with smoothing
plot(poly_trim)

# Maximum smooth is at 50 (so smooth)
poly_smooth <- smoothr::smooth(poly_trim, method = "ksmooth", smoothness = 20)
plot(poly_smooth)


poly_smooth %>%
  ggplot() +
  stars::geom_stars(
    data = dem_filter %>% 
      raster::trim() %>% 
      stars::st_as_stars()
  ) +
  geom_sf(data = poly_trim, color = "grey50", fill = "transparent") +
  geom_sf(fill = "transparent", size = 1) +
  scale_fill_distiller(palette = "Spectral")


# Write polygon shape to file

# Project in standard coordinate system
poly_smooth <- poly_smooth %>%
  # No attribute fields are necessary at this point
  dplyr::summarize() %>%
  sf::st_transform(crs = spatial_ref$epsg)

if (!dir.exists(file.path(path_out, "delineation"))) {
  dir.create(file.path(path_out, "delineation"))
} 
sf::st_write(
  poly_smooth, 
  file.path(path_out, "delineation", paste0(settings$site, "_delin.shp")), 
  delete_layer = TRUE
)

