
rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

control <- list(
  glcm = c("rgb_1", "rgb_2", "rgb_3", "nir")
)

# References

# Haralick, R. M., Shanmugam, K., & Dinstein, I. (1973). Textural Features for 
# Image Classification. IEEE Transactions on Systems, Man, and Cybernetics, 
# SMC-3(6), 610–621. https://doi.org/10.1109/TSMC.1973.4309314
# 
# van der Werff, H. M. A., & van der Meer, F. D. (2008). Shape-based 
# classification of spectrally identical objects. ISPRS Journal of 
# Photogrammetry and Remote Sensing, 63(2), 251–258. 
# https://doi.org/10.1016/j.isprsjprs.2007.09.007


library(progress)
library(tidyverse)

raster_rescale <- function(x, ...) {
  values <- raster::values(x)
  new_values <- scales::rescale(values, ...)
  raster::setValues(x, new_values)
}

st_width <- function(x) {
  bbox <- sf::st_bbox(x)
  unname(bbox$xmax - bbox$xmin)
}

st_height <- function(x) {
  bbox <- sf::st_bbox(x)
  unname(bbox$ymax - bbox$ymin)
}

features_basic <- function(segm, feat) {
  
  basic <- EBImage::computeFeatures.basic(segm, feat, basic.quantiles = NA)
  
  basic %>%
    tibble::as_tibble(rownames = "segment") %>%
    dplyr::select(1:3) %>%
    dplyr::rename_with(~ stringr::str_replace(.x, "\\.", "_")) %>%
    dplyr::rename_with(~ stringr::str_remove(.x, "^b_")) %>%
    dplyr::mutate(segment = as.integer(segment))
}

features_glcm <- function(segm, feat) {
  
  # From the original reference, Haralick et al. 1973
  # 1. Angular second moment (asm)
  # 2. Contrast (con)
  # 3. Correlation (cor)
  # 4. Variance (var)
  # 5. Inverse difference moment (idm)
  # 6. Sum average (sav)
  # 7. Sum variance (sva)
  # 8. Sum entropy (sen)
  # 9. Entropy (ent)
  # 10. Difference variance (dva)
  # 11. Difference entropy (den)
  # 12. Information measure of correlation 1 (f12)
  # 13. Information measure of correlation 2 (f13)
  
  glcm <- EBImage::computeFeatures.haralick(segm, feat, haralick.scales = 1)
  
  glcm %>% 
    tibble::as_tibble(rownames = "segment") %>% 
    dplyr::rename_with(~ stringr::str_remove(.x, ".s1")) %>%
    dplyr::rename_with(~ stringr::str_replace(.x, "\\.", "_")) %>%
    dplyr::rename_with(~ stringr::str_remove(.x, "^h_")) %>%
    dplyr::mutate(segment = as.integer(segment)) %>%
    dplyr::filter(segment %in% unique(segm_mat))
}

progress_info <- function(len) {
  
  progress_bar$new(
    total = len, 
    format = paste0(
      "[:spin] Completed: :current | :percent  ", 
      "Elapsed: :elapsed  Remaining: :eta"
    ),
    clear = FALSE
  )
}


### Inputs

# Set working directory
wd <- file.path(
  "/Users/Graham/Desktop/DATA/Flux", settings$site, "analysis", "spatial"
)

paths <- list(
  spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  # Path to image segments
  segm = file.path(wd, "04_segmentation", "segments_rgbn.tif"),
  # Path to raster data
  feat = file.path(wd, "03_feature_extraction"),
  # Path to output file
  out = file.path(wd, "05_segment_features", "segment_features.csv")
)

# Load reference files
source(paths$spatial_ref)


### Load input data ============================================================

# Load image segments
segm <- paths$segm %>% 
  stars::read_stars() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::select(segment = 1) %>%
  dplyr::mutate(segment = as.integer(segment))

# Import raster data
feat_rst <- paths$feat %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE) %>%
  raster::stack()

# Read raster data to tidy format
# - this is a hacky way of doing it but c.stars simply doesn't work (why??) 
feat_st <- paths$feat %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE) %>%
  purrr::map(stars::read_stars) %>% 
  purrr::map(tibble::as_tibble, center = FALSE) %>% 
  purrr::map_if(
    ~ "band" %in% names(.), 
    ~ tidyr::pivot_wider(
      .x, names_from = band, values_from = 4, names_glue = "{.value}_{band}"
    )
  ) %>% 
  purrr::reduce(dplyr::full_join, by = c("x", "y")) %>%
  dplyr::rename_with(~ stringr::str_remove(.x, ".tif")) %>%
  stars::st_as_stars() %>%
  sf::st_set_crs(spatial_ref$epsg)

### Assign features to segments ================================================

# Convert segments to polygon format 
segm_poly <- sf::st_as_sf(segm, as_points = FALSE, merge = TRUE)

# Calculate segment shape features
# - via van der Werff & van der Meer 2008
data_shape <- segm_poly %>%
  dplyr::mutate(
    height = purrr::map_dbl(geometry, st_height),
    width = purrr::map_dbl(geometry, st_width),
    area = sf::st_area(geometry),
    perimeter = lwgeom::st_perimeter(geometry),
    convex_perim = geometry %>% sf::st_convex_hull() %>% lwgeom::st_perimeter(),
    compactness = 4 * pi * area / perimeter^2,
    roundness = 4 * pi * area / convex_perim^2,
    convexity = convex_perim / perimeter
  ) %>%
  sf::st_drop_geometry() %>%
  tibble::as_tibble() %>%
  dplyr::select(-convex_perim) %>%
  dplyr::mutate(dplyr::across(where(~ inherits(.x, "units")), as.numeric))

# Calculate segment basic & texture features
# - might be slow, so set up loop to check progress
n_feat <- length(feat_st)
data_feat_list <- vctrs::vec_init(list(), n_feat)
segm_mat <- segm %>% purrr::pluck(1) %>% t()
p <- progress_info(n_feat)

for (i in 1:n_feat) {
  
  name_temp <- names(feat_st[i])
  feat_temp <-  raster::as.matrix(feat_rst[[i]])
  
  data_temp <- features_basic(segm_mat, feat_temp)
  
  if (name_temp %in% control$glcm) {
    # GLCM input needs to be in [0, 1] range
    glcm_temp <- features_glcm(segm_mat, scales::rescale(feat_temp))
    
    data_temp <- dplyr::full_join(data_temp, glcm_temp, by = "segment")
  }
  
  data_feat_list[[i]] <- dplyr::rename_with(
    data_temp, ~ stringr::str_c(name_temp, "_", .x), -segment
  )
  
  p$tick()
}

# Bind together all feature data
data_feat <- data_feat_list %>%
  purrr::reduce(dplyr::full_join, by = "segment") %>%
  dplyr::full_join(data_shape, by = "segment") %>%
  dplyr::arrange(segment)

# Write to file
readr::write_csv(data_feat, paths$out)


# Extract mean & variance for all features in each image segment
# - takes a couple minutes, be patient
# feature_data <- segments_poly %>%
#   dplyr::select(segment, geometry) %>%
#   sf::st_join(features_st, ., what = "inner", as_points = TRUE) %>%
#   sf::st_drop_geometry() %>%
#   tibble::as_tibble() %>%
#   dplyr::select(-x, -y) %>%
#   dplyr::group_by(segment) %>%
#   # Using scoped summarize because across is WAY too slow (dplyr 1.0.0)
#   dplyr::summarize_all(list(mean = base::mean, sd = stats::sd)) %>%
#   # Bind polygon-level metrics & geometry
#   dplyr::left_join(segments_poly, by = "segment") %>%
#   sf::st_as_sf() %>%
#   sf::st_join(., train_points, what = "inner", as_points = TRUE) %>%
#   sf::st_drop_geometry() %>%
#   tibble::as_tibble()

