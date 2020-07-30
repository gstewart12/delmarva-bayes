
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLN", # three letter site code
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

# References

# Balázs, B., Bíró, T., Dyke, G., Singh, S. K., & Szabó, S. (2018). Extracting 
# water-related features using reflectance data and principal component analysis 
# of Landsat images. Hydrological Sciences Journal, 63(2), 269–284. 
# https://doi.org/10.1080/02626667.2018.1425802
# 
# Dabrowska-Zielinska, K., Gruszczynska, M., Lewinski, S., Hoscilo, A., & 
# Bojanowski, J. (2009). Application of remote and in situ information to the 
# management of wetlands in Poland. Journal of Environmental Management, 90(7), 
# 2261–2269. https://doi.org/10.1016/j.jenvman.2008.02.009
# 
# Lefebvre, G., Davranche, A., Willm, L., Campagna, J., Redmond, L., Merle, C., 
# et al. (2019). Introducing WIW for Detecting the Presence of Water in Wetlands 
# with Landsat and Sentinel Satellites. Remote Sensing, 11(19), 2210. 
# https://doi.org/10.3390/rs11192210
# 
# Luo, S., Wang, C., Xi, X., Pan, F., Qian, M., Peng, D., et al. (2017). 
# Retrieving aboveground biomass of wetland Phragmites australis (common reed) 
# using a combination of airborne discrete-return LiDAR and hyperspectral data. 
# International Journal of Applied Earth Observation and Geoinformation, 58, 
# 107–117. https://doi.org/10.1016/j.jag.2017.01.016
# 
# Maguigan, M., Rodgers, J., Dash, P., & Meng, Q. (2016). Assessing Net Primary 
# Production in Montane Wetlands from Proximal, Airborne, and Satellite Remote 
# Sensing. Advances in Remote Sensing, 5(2), 118–130. 
# https://doi.org/10.4236/ars.2016.52010
# 
# Millard, K., & Richardson, M. (2015). On the Importance of Training Data 
# Sample Selection in Random Forest Image Classification: A Case Study in 
# Peatland Ecosystem Mapping. Remote Sensing, 7(7), 8489–8515. 
# https://doi.org/10.3390/rs70708489
# 
# Taddeo, S., Dronova, I., & Depsky, N. (2019). Spectral vegetation indices of 
# wetland greenness: Responses to vegetation structure, composition, and spatial 
# distribution. Remote Sensing of Environment, 234, 111467. 
# https://doi.org/10.1016/j.rse.2019.111467


library(tidyverse)

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

calc_gemi <- function(r, n) {
  ngemi <- 2 * (n^2 - r^2) + 1.5 * n + 0.5 * r
  nr <- n + r
  ngemi / (nr + 0.5) * (1 - 0.25 * ngemi / (nr + 0.5)) - (r - 0.125) / (1 - r)
}

fill_raster_clumps <- function(x, size_thr = 5, fill_value = 0, na = NA) {
  
  x[raster::values(x) %in% na] <- NA
  
  n_na <- length(x[is.na(x)])
  
  if (n_na == 0) {
    return(x)
  }
  
  focal_fill <- raster::focal(
    x, w = matrix(1, 3, 3), fun = median, na.rm = TRUE, pad = TRUE, 
    NAonly = TRUE
  )
  
  focal_fill_ind <- x %>% 
    # Reclassify raster as NA or non-NA
    raster::reclassify(cbind(-Inf, Inf, 0)) %>% 
    raster::reclassify(cbind(NA, 1)) %>% 
    # "Clump" connected cells
    raster::clump() %>% 
    raster::values() %>% 
    tibble::enframe(name = "cell", value = "clump") %>% 
    tidyr::drop_na() %>% 
    # Number of cells in each clump 
    dplyr::group_by(clump) %>% 
    dplyr::mutate(size = dplyr::n()) %>% 
    dplyr::ungroup() %>%
    # Identify cells within clumps smaller than size threshold
    dplyr::filter(size < size_thr) %>%
    dplyr::pull(cell)
  
  out <- x
  # Fill 'eligable' cells with focal median interpolation
  out[focal_fill_ind] <- focal_fill[focal_fill_ind]
  # Fill the rest with a value
  out[is.na(out)] <- fill_value
  out
}

# Set working directory
wd <- file.path(
  "/Users/Graham/Desktop/DATA/Flux", settings$site, "analysis", "spatial"
)

paths <- list(
  spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  # Set site delineation path
  delin = file.path(wd, "02_site_delineation", "delineation"),
  # Set rgbn image path
  rgbn =  file.path(wd, "01_clip_imagery", "naip_2018.tif"),
  # Set DEM path
  dem = file.path(wd, "01_clip_imagery", "dem_2007.tif"),
  # Set LiDAR path
  lidar = file.path(wd, "01_clip_imagery", "lidar_2014.las"),
  # Set output path
  out = file.path(wd, "03_feature_extraction")
)

# Load reference files
source(paths$spatial_ref)

# Read site delineation
delin <- sf::read_sf(paths$delin)
site_bbox <- sf::st_bbox(delin)


### Spectral data ==============================================================

# Color spaces
color <- list()

# Read RGBN data
rgbn_st <- paths$rgbn %>% 
  stars::read_stars() %>% 
  sf::st_set_crs(spatial_ref$epsg) %>%
  stars::st_set_dimensions("band", values = c("r", "g", "b", "n")) 

rgbn <- rgbn_st %>%
  tibble::as_tibble(center = FALSE) %>% 
  tidyr::pivot_wider(names_from = band, values_from = 4)

# Near-infrared (extract this from RGBN image)
color$nir <- dplyr::select(rgbn, x, y, nir = n)

# RGB
color$rgb <- dplyr::select(rgbn, -n)

# HSV
color$hsv <- color$rgb %>%
  dplyr::bind_cols(tibble::as_tibble(t(grDevices::rgb2hsv(.$r, .$g, .$b)))) %>%
  # Rescale to "spectral range"
  dplyr::mutate(dplyr::across(c(h, s, v), ~ .x * 255)) %>%
  dplyr::select(-c(r:b))

# XYZ
# color$xyz <- color$rgb %>% 
#   dplyr::select(r, g, b) %>% 
#   grDevices::convertColor(
#     from = "sRGB", to = "XYZ", scale.in = 255, scale.out = 255
#   ) %>%
#   tibble::as_tibble(.name_repair = "minimal") %>%
#   rlang::set_names(c("X", "Y", "Z")) %>%
#   dplyr::bind_cols(dplyr::select(color$rgb, x, y), .)

# C1C2C3 - invariant to shadowing effects
color$c1c2c3 <- color$rgb %>% 
  dplyr::rowwise() %>%
  dplyr::mutate(
    # Extended gradual C1C2C3 (Golchin et al. 2013) - max = 1
    #c1 = atan(r / sqrt(g^2 + b^2)),
    #c2 = atan(g / sqrt(r^2 + b^2)),
    #c3 = atan(b / sqrt(r^2 + g^2))
    # Original C1C2C3 (Sun & Li 2010) - max = atan(255)
    c1 = atan(r / max(g, b)),
    c2 = atan(g / max(r, b)),
    c3 = atan(b / max(r, g))
  ) %>%
  dplyr::ungroup() %>%
  # Rescale to "spectral range"
  dplyr::mutate(dplyr::across(c(c1, c2, c3), ~ .x / atan(255) * 255)) %>%
  dplyr::select(-c(r:b))

# l1l2l3 - invariant to highlighting effects
color$l1l2l3 <- color$rgb %>% 
  dplyr::mutate(
    l1 = (r - g)^2 / ((r - g)^2 + (r - b)^2 + (g - b)^2),
    l2 = (r - b)^2 / ((r - g)^2 + (r - b)^2 + (g - b)^2),
    l3 = (g - b)^2 / ((r - g)^2 + (r - b)^2 + (g - b)^2)
  ) %>%
  # If r == b == g then l1l2l3 will be NaN (fill these values)
  dplyr::mutate(dplyr::across(l1:l3, ~ tidyr::replace_na(.x, 1 / 3))) %>%
  # Rescale to "spectral range"
  dplyr::mutate(dplyr::across(c(l1, l2, l3), ~ .x * 255)) %>%
  dplyr::select(-c(r:b))

# Add spatial reference to color matrices
color_st <- color %>%
  purrr::map(stars::st_as_stars) %>%
  purrr::map(sf::st_set_crs, spatial_ref$epsg) %>%
  purrr::map(~ stars::st_redimension(.x, along = list(band = names(.x)))) %>%
  # Crop to site area limits
  purrr::map(sf::st_crop, site_bbox, as_points = FALSE)

# Write color spaces to file
color_out <- file.path(paths$out, "color")
unlink(paste0(color_out, "/*")) # clear folder contents (CANNOT BE UNDONE)
purrr::imap(
  color_st, ~ stars::write_stars(.x, file.path(color_out, paste0(.y, ".tif")))
)


# Spectral indices
# - chosen based on inclusion in wetland remote sensing studies
# - specifically, Dabrowska-Zielinska et al. 2009, Maguigan et al. 2016, 
#   Balazs et al. 2018, Lefebvre et al. 2019, Taddeo et al. 2019

# Shadow-eliminated Vegetation Index (Jiang et al. 2019)
# - need to optimize a parameter before calculation 
x_dim <- rgbn_st %>% stars::st_get_dimension_values(which = 1) %>% length()
f_delta <- seq(0, 1 * 255, by = 0.001 * 255) %>%
  # Sample from the first 10 rows 
  purrr::map(~ dplyr::mutate(
    rgbn[1:(10 * 300), ], rvi = n / r, svi = 1 / r, sevi = rvi + .x * svi
  )) %>% 
  purrr::map(~ dplyr::summarize(
    .x, r1 = cor(sevi, rvi), r2 = cor(sevi, svi), f = abs(r1 - r2)
  )) %>%
  purrr::map_dbl(dplyr::pull, f) %>%
  rlang::set_names(seq(0, 1 * 255, by = 0.001 * 255)) %>% 
  tibble::enframe() %>% 
  dplyr::arrange(value) %>%
  dplyr::transmute(name = as.numeric(name)) %>%
  purrr::pluck("name", 1)

index <- dplyr::transmute(
  rgbn,
  
  # VEGETATION INDICES
  
  # Simple Ratio (0, 255)
  sr = r/n,
  # Ratio Vegetation Index - Richardson and Wiegand 1977 (0, 255) 
  rvi = n/r,
  # Normalized Difference Vegetation Index - Rouse et al. 1974 (-1, 1)
  ndvi = (n - r) / (n + r),
  # Renormalized Difference Vegetation Index ~(-16, 16) 
  rdvi = (n - r) / sqrt(n + r),
  # Normalized Ratio Vegetation Index - Baret & Guyot 1991 (-1, 1)
  nrvi = (r/n - 1) / (r/n + 1),
  # Modified Simple Ratio - Chen & Cihlar 1996 (-1, 1)
  msr = (n/r - 1) / (sqrt(n/r) + 1),
  # Thiam's Transformed Vegetation Index - Thiam 1997 (0, 1)
  ttvi = sqrt(abs(ndvi + 0.5)),
  # Green Normalized Difference Vegetation Index (-1, 1)
  gndvi = (n - g) / (n + g),
  # Wide Dynamic Range Vegetation Index - Gitelson 2004 (-1, 1) 
  wdrvi = (0.1 * n - r) / (0.1 * n + r),
  # Atmospherically Resistent Vegetation Index - Kaufman & Tanre 1992 (-1, 1)
  arvi = (n - (2 * r - b)) / (n + (2 * r - b)),
  # Difference Vegetation Index (-255, 255)
  dvi = 1 * n - r,
  # Soil Adjusted Vegetation Index - Huete 1988 (-1.5, 1.5)
  savi = 1.5 * (n - r) / (n + r + 0.5),
  # Modified Soil Adjusted Vegetation Index
  # - used in Luo et al. 2017
  msavi = n + 0.5 - (0.5 * sqrt((2 * n + 1)^2 - 8 * (n - (2 * r)))),
  # Corrected Transformed Vegetation Index - Perry and Lautenschlager 1984
  # - used by Balazs et al. 2018
  # - as (ndvi + 0.5) / abs(ndvi + 0.5) * sqrt(abs(ndvi + 0.5))
  ctvi = (n + 0.5) / sqrt(abs(n + 0.5)),
  # Visible Atmospherically Resistant Index - Gitelson et al. 2002
  vari = (g - r) * (g + r - b),
  # Green Chromatic Coordinate (0, 1)
  gcc = g / (g + r + b),
  # Global Environment Monitoring Index - Pinty & Verstraete 1992
  gemi = calc_gemi(r, n),
  # Shadow Eliminated Vegetation Index - Jiang et al. 2019
  sevi = n/r + f_delta * 1/r,
  
  # WATER INDICES - via Malinowski et al. 2015
  # Normalized Difference Water Index - McFeeters 1996 (-1, 1)
  ndwi = (g - n) / (g + n),
  # Difference between Vegetation and Water - Gond et al. 2004
  dvw = ndvi - ndwi,
  # Index of Free Water - Adell & Puech 2003 (-255, 255)
  ifw = n - g,
  # Water Index (0, 65025)
  wi = n^2 / b,
  # Water Impoundment Index - Caillaud et al. 1987 (0, 65025)
  wii = n^2 / r,
  # Normalized Difference Turbidity Index - Lacaux et al. 2007 (-1, 1)
  ndti = (r - g) / (r + g)
)

# Add spatial reference to spectral indices
index_st <- index %>% 
  purrr::map(tibble::as_tibble) %>% 
  purrr::imap(~ rlang::set_names(.x, .y)) %>%
  purrr::map(~ dplyr::bind_cols(dplyr::select(rgbn, x, y), .x)) %>%
  purrr::map(stars::st_as_stars) %>%
  purrr::map(sf::st_set_crs, spatial_ref$epsg) %>%
  purrr::map(~ stars::st_redimension(.x, along = list(band = names(.x)))) %>%
  # Crop to site area limits
  purrr::map(sf::st_crop, site_bbox, as_points = FALSE)

# Write spectral indices to file
index_out <- file.path(paths$out, "index")
unlink(paste0(index_out, "/*")) # clear folder contents (CANNOT BE UNDONE)
purrr::imap(
  index_st, ~ stars::write_stars(.x, file.path(index_out, paste0(.y, ".tif")))
)


### LiDAR data =================================================================

# Read DEM raster (already fully processed)
dem <- paths$dem %>% 
  stars::read_stars() %>% 
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::select(dem = 1)

# Topographic metrics are not very useful here due to low gradient

# DEM will be a template to resample grid calculations
grid <- as(dem, "Raster")

# Read raw LiDAR points
las <- lidR::readLAS(paths$lidar)

# Compute digital terrain model (DTM)
dtm <- lidR::grid_terrain(las, algorithm = lidR::knnidw(), res = grid)

# Compute height above ground (HAG)
hag <- lidR::normalize_height(las, algorithm = lidR::knnidw())
hag$Z[hag$Z < 0] <- 0 # fix negative heights

# Helper function for gridded HAG metrics

safe_max <- function(..., na.rm = FALSE) {
  max <- purrr::quietly(max)(..., na.rm = na.rm)$result
  dplyr::na_if(max, -Inf)
} 
safe_min <- function(..., na.rm = FALSE) {
  min <- purrr::quietly(min)(..., na.rm = na.rm)$result
  dplyr::na_if(min, Inf)
} 

metric_fns <- function(z, i, rn, class, dz = 0.01) {
  
  # Most of these are explained in Millard & Richardson 2015 or Luo et al. 2017
  
  first <- rn == 1
  ground <- class == 2
  veg <- class != 2
  n <- length(z)
  i_sum <- sum(i)
  
  z_metrics <- list(
    z_max = max(z), 
    z_min = min(z), 
    z_mean = mean(z),
    z_sd = sd(z),
    z_cv = sd(z) / mean(z),
    z_vmin = safe_min(z[veg]),
    z_vmean = mean(z[veg]),
    z_vsd = sd(z[veg]),
    z_vcv = sd(z[veg]) / mean(z[veg]),
    z_1mean = mean(z[first]),
    zi_mean = mean(z * i)
    # Entropy and VCI calibrated for very low vegetation
   # z_entropy = if (all(ground)) 0 else lidR::entropy(z, by = dz),
    #vci = lidR::VCI(z, by = dz, zmax = 30)
  )
  
  i_metrics <- list(
    i_max = max(i),
    i_min = min(i),
    i_mean = mean(i),
    i_sd = sd(i),
    i_cv = sd(i) / mean(i),
    i_1mean = mean(i[first]),
    i_gmean = mean(i[ground]),
    i_1p = sum(i[first]) / i_sum,
    # Percent of canopy intensity, a description of canopy density
    i_vp = sum(i[veg]) / i_sum
  )
  
  rn_metrics <- list(
    rn_max = max(rn),
    rn_mean = mean(rn),
    rn_1p = sum(first) / n,
    # Percent of canopy returns, a description of canopy cover
    rn_vp = sum(veg) / n
  )
  
  c(z_metrics, i_metrics, rn_metrics)
} 

.metric_fns <- ~ metric_fns(Z, Intensity, ReturnNumber, Classification)

# Calculate LiDAR metrics in raster format
# - for some reason inputting grid as res doesn't work 
lidar_metrics <- lidR::grid_metrics(hag, func = .metric_fns, res = grid)
names(lidar_metrics)
# Resample to site grid/coordinate system
#lidar_metrics_res <- raster::resample(lidar_metrics, grid, method = "ngb")

# Fill missing values
# - general assumption is that NA means no returns were measured
# - most likely due to signal degradation through water
# - therefore 0 returns = 0 intensity, 0 height above ground -> fill_value = 0
# - exception: NAs in rn_1p, rn_max, rn_mean filled with 1 
#   - assumes that if any pulses returned, would have been ground and only 1
# - exception: NAs in i_1p filled with 1
#   - assumes all of the intensity (i.e. 0) would have been in first return
# - exception: NAs in *_cv filled with 1

# As of lidR 3.0.1 need to manually set NA values
na_metrics <- lidar_metrics[["rn_max"]] == 0

# Fill value key
fill_value_key <- lidar_metrics %>%
  names() %>%
  rlang::set_names(~ .) %>%
  purrr::map(~ 0) %>%
  purrr::modify_at(
    c("rn_max", "rn_mean", "rn_1p", "i_1p", "z_cv", "z_vcv", "i_cv"), ~ 1
  )

# Perform gap filling
lidar <- lidar_metrics
lidar[na_metrics] <- NA
lidar <- lidar %>% 
  raster::unstack() %>%
  rlang::set_names(names(lidar)) %>%
  purrr::map2(fill_value_key, ~ fill_raster_clumps(.x, fill_value = .y))

# Standard format
lidar_st <- lidar %>%
  purrr::map(stars::st_as_stars) %>%
  purrr::map(sf::st_set_crs, spatial_ref$epsg) %>%
  # Add DEM to write with the rest of the LiDAR metrics
  purrr::prepend(list(dem = dem)) %>%
  purrr::map(~ stars::st_redimension(.x, along = list(band = names(.x)))) %>%
  # Crop to site area limits
  purrr::map(sf::st_crop, site_bbox, as_points = FALSE)

# Write LiDAR metrics to file
lidar_out <- file.path(paths$out, "lidar")
unlink(paste0(lidar_out, "/*")) # clear folder contents (CANNOT BE UNDONE)
purrr::imap(
  lidar_st, ~ stars::write_stars(.x, file.path(lidar_out, paste0(.y, ".tif")))
)

# Finished

