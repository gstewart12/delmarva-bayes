
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)


library(tidyverse)

# Outline:
# 1. Load imagery
# 2. Transform coordinate system
# 3. Crop to site area
# 4. Resample to site grid
# 5. Apply imagery-specific corrections

wd <- "/Users/Graham/Desktop/DATA/Spatial"

paths <- list(
  spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  aoi = file.path(dirname(wd), "Flux", settings$site, "site_info", "aoi"),
  point = file.path(wd, "survey", "flux_pts"),
  naip = file.path(wd, "NAIP", "m_3907558_se_18_060_20180810.tif"),
  dem = file.path(wd, "DEM", "2007_1m_DEM_filled.tif"),
  lidar = file.path(
    wd, "LiDAR", "de_md2014_usgs_cmgp_sandy", "de_md2014_usgs_cmgp_sandy.las"
  ),
  out = file.path(
    dirname(wd), "Flux", settings$site, "analysis", "spatial", "01_clip_imagery"
  )
)

# Load reference files
source(paths$spatial_ref)

plot_stars <- function(data, band, ..., rgb = FALSE) {
  
  band <- rlang::enquo(band)
  
  if (rgb) {
    ggplot2::ggplot() + 
      geom_stars_rgb(data = data, ...) +
      ggplot2::scale_fill_identity() +
      ggplot2::coord_equal() +
      ggplot2::theme_void()
  } else {
    if (!rlang::quo_is_missing(band)) data <- dplyr::select(data, !!band)
    
    ptype <- data %>% 
      dplyr::pull(1) %>% 
      vctrs::vec_ptype()
    
    plot <- ggplot2::ggplot() + 
      stars::geom_stars(data = data) +
      ggplot2::coord_equal() +
      ggplot2::theme_void()
    
    if (is.factor(ptype)) {
      plot + ggplot2::scale_fill_viridis_d(option = "A")
    } else {
      plot + ggplot2::scale_fill_viridis_c(option = "A")
    }
  }
}

geom_stars_rgb <- function(data, r = 1, g = 2, b = 3, max_value = 255, ...) {
  
  if (length(stars::st_dimensions(data)) > 2) {
    crs <- sf::st_crs(data)
    data <- data %>% 
      tibble::as_tibble(center = FALSE) %>% 
      tidyr::pivot_wider(names_from = 3, values_from = 4) %>%
      stars::st_as_stars() %>%
      sf::st_set_crs(crs)
  }
  
  data <- data %>%
    dplyr::select(
      r = dplyr::all_of(r), g = dplyr::all_of(g), b = dplyr::all_of(b)
    ) %>%
    dplyr::mutate(rgb = grDevices::rgb(r, g, b, maxColorValue = max_value))
  stars::geom_stars(data = data, ggplot2::aes(x = x, y = y, fill = rgb), ...)
}


### AOI initialization =========================================================

# Read AOI
aoi <- sf::read_sf(paths$aoi)
aoi <- sf::st_set_crs(aoi, spatial_ref$epsg)

# Read points of interest
point <- sf::read_sf(paths$point)
tower <- point %>% 
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::filter(site == settings$site, type == "tower") %>%
  sf::st_geometry()

# Ensure that grid template will cover entire AOI
# - square buffer of [length] around tower should contain AOI regardless of 
#   exact tower location
length <- aoi %>% 
  sf::st_geometry() %>% 
  sf::st_area() %>% 
  as.numeric() %>% 
  sqrt() %>% 
  ceiling()
  # Center cells on tower (so tower point is in middle of cell, not on vertex)
  #magrittr::add(0.5)

# Construct grid template
grid <- tower %>% 
  # Create large box centered on tower
  sf::st_buffer(length) %>% 
  sf::st_bbox() %>% 
  # Rasterize
  stars::st_as_stars(nx = length * 2, ny = length * 2) %>%
  # Trim area, allowing extra room for focal stats to be calculated
  sf::st_crop(sf::st_bbox(sf::st_buffer(aoi, 15)), as_points = FALSE)

# Make sure tower is where expected on grid
grid %>% 
  sf::st_as_sf() %>%
  sf::st_crop(sf::st_buffer(tower, 5)) %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf() +
  ggplot2::geom_sf(data = tower, size = 2)



### NAIP aerial imagery ========================================================

# Read imagery
naip <- raster::stack(paths$naip)

# Crop and snap to site grid
naip_crop <- raster::resample(naip, as(grid, "Raster"), method = "bilinear")

# Can't use st_warp for this because no support for bilinear interpolation
# naip_crop <- naip %>%
#   stars::st_as_stars() %>%
#   stars::st_warp(grid)
  
# Make prettier band names
naip_bands <- naip_crop %>%
  stars::st_as_stars() %>%
  stars::st_set_dimensions(
    "band", values = c("red", "green", "blue", "nir")
  ) %>%
  dplyr::select(value = 1)

# Check result
plot_stars(naip_bands, rgb = TRUE)

# Save site-level imagery
stars::write_stars(naip_bands, file.path(paths$out, "naip_2018.tif"))


### Digital elevation model (DEM) ==============================================

# Read DEM
dem <- raster::raster(paths$dem)

# Invert values
# dem_min <- raster::minValue(dem)
# dem_max <- raster::maxValue(dem)
# dem_invert <- (dem - dem_max) * -1 + dem_min

# Crop and snap to site grid
dem_crop <- dem %>%
  raster::resample(as(grid, "Raster"), method = "bilinear") %>%
  stars::st_as_stars()

# Crop and snap to site grid
# dem_crop <- dem_invert %>%
#   stars::st_as_stars() %>%
#   stars::st_warp(grid)

# Check result
plot_stars(dem_crop)

# Save site-level DEM
stars::write_stars(dem_crop, file.path(paths$out, "dem_2007.tif"))


### LiDAR points ===============================================================

# Read LiDAR points
lidar <- lidR::readLAS(paths$lidar)

# Crop and snap to site grid
lidar_proj <- sp::spTransform(lidar, spatial_ref$proj4)

lidar_bbox <- grid %>% 
  as("Raster") %>%
  raster::extent() %>%
  # Expand area slightly since not rasterizing yet
  magrittr::add(2)
  
lidar_crop <- lidR::clip_roi(lidar_proj, lidar_bbox)

# Save site-level LiDAR data
lidR::writeLAS(lidar_crop, file.path(paths$out, "lidar_2014.las"))


