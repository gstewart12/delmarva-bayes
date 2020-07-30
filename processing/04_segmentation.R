
rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

# References

# De Luca, G., N. Silva, J. M., Cerasoli, S., Araújo, J., Campos, J., Di Fazio, 
# S., & Modica, G. (2019). Object-Based Land Cover Classification of Cork Oak 
# Woodlands using UAV Imagery and Orfeo ToolBox. Remote Sensing, 11(10), 1238. 
# https://doi.org/10.3390/rs11101238
# 
# Deborah, H., Richard, N., & Hardeberg, J. Y. (2015). A Comprehensive 
# Evaluation of Spectral Distance Functions and Metrics for Hyperspectral Image 
# Processing. IEEE Journal of Selected Topics in Applied Earth Observations and 
# Remote Sensing, 8(6), 3224–3234. https://doi.org/10.1109/JSTARS.2015.2403257
# 
# Michel, J., Youssefi, D., & Grizonnet, M. (2015). Stable Mean-Shift Algorithm 
# and Its Application to the Segmentation of Arbitrarily Large Remote Sensing 
# Images. IEEE Transactions on Geoscience and Remote Sensing, 53(2), 952–964. 
# https://doi.org/10.1109/TGRS.2014.2330857
# 
# Ming, D., Li, J., Wang, J., & Zhang, M. (2015). Scale parameter selection by 
# spatial statistics for GeOBIA: Using mean-shift based multi-scale segmentation 
# as an example. ISPRS Journal of Photogrammetry and Remote Sensing, 106, 28–41. 
# https://doi.org/10.1016/j.isprsjprs.2015.04.010


control <- list(
  # Files to extract imagery from
  # - USE HSV INSTEAD; BETTER AT DETECTING OBJECTS? (Rasouli & Tsotsos 2017)
  img_sources = c("color/rgb", "color/nir"),
  # Raster layers used in segmentation (really only used for naming)
  bands = c("r", "g", "b", "nir"),
  # Mean shift & segmentation parameters (OTB default)
  # - ranger: spectral "size" or compactness of classes (15)
  # - spatialr: spatial size of the patches making up a class (5)
  # - minsize: minimum acceptable spatial size of patches in a class (0)
  maxiter = 100,
  minsize = 5
  #mss_pars = list(ranger = 15, spatialr = 5),
  #lsms_pars = list(ranger = 15, spatialr = 5)
)

#library(progress)
library(sf)
library(tidyverse)


raster_rescale <- function(x, ...) {
  values <- raster::values(x)
  new_values <- scales::rescale(values, ...)
  raster::setValues(x, new_values)
}

plot_stars <- function(data, rgb = FALSE, ...) {
  
  if (rgb) {
    ggplot2::ggplot() + 
      geom_stars_rgb(data = data, ...) +
      ggplot2::scale_fill_identity() +
      ggplot2::coord_equal() +
      ggplot2::theme_void()
  } else {
    ggplot2::ggplot() + 
      stars::geom_stars(data = data) +
      ggplot2::scale_fill_viridis_c(option = "A") +
      ggplot2::coord_equal() +
      ggplot2::theme_void()
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

lv <- function(x, w = 3) {
  #raster::focal(x, w = matrix(1, w, w), fun = stats::sd) 
  # Pretty much the same as focal sd, but way faster
  mu <- raster::focal(x, w = matrix(1 / w^2, w, w), fun = sum)
  sqrt(raster::focal((x - mu)^2, w = matrix(1 / w^2, w, w), fun = sum))
}

alv <- function(x, w = 3) {
  variance <- lv(x, w = w)
  raster::cellStats(variance, mean, na.rm = TRUE)
}

# Set working directory
wd <- file.path(
  "/Users/Graham/Desktop/DATA/Flux", settings$site, "analysis", "spatial"
)

paths <- list(
  spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  md = "/Users/Graham/Desktop/DATA/Flux/tools/reference/site_metadata.R",
  # Path to site delineation
  delin = file.path(wd, "02_site_delineation", "delineation"),
  # Path to raster data used for image segmentation
  feat = file.path(wd, "03_feature_extraction"),
  # Output path
  out = file.path(wd, "04_segmentation")
)

# Temporary file paths
temp_files <- list(
  img_stack = file.path(tempdir(), "img_stack.tif"),
  filter_range = file.path(tempdir(), "otb_filter_range.tif"),
  filter_spatial = file.path(tempdir(), "otb_filter_spatial.tif"),
  segm_init = file.path(tempdir(), "otb_segm_init.tif"),
  segm_merge = file.path(tempdir(), "otb_segm_merge.tif")
)

# File paths for all outputs
files_out <- list(
  params = file.path(paths$out, "segm_params.txt"),
  rst_init = file.path(paths$out, "segments_init_rgbn.tif"),
  rst = file.path(paths$out, "segments_rgbn.tif"),
  poly = file.path(paths$out, "segments", "segments_rgbn.shp")
)

# Load reference files
source(paths$spatial_ref)
source(paths$md)
md <- purrr::pluck(site_metadata, settings$site)

### Load input data

# Set tower location
tower <- sf::st_sfc(sf::st_point(c(md$x_utm, md$y_utm)), crs = spatial_ref$epsg)

# Read AOI
delin <- sf::read_sf(paths$delin)

# Read imagery, format, and save to temp dir
# - layers used for segmentation should provide information at the desired scale 
#   of individual segments
img_stack <- paths$feat %>%
  list.files(pattern = ".tif$", full.name = TRUE, recursive = TRUE) %>%
  stringr::str_subset(
    stringr::str_c(control$img_sources, ".", collapse = "|")
  ) %>%
  magrittr::extract(order(control$img_sources)) %>%
  purrr::map(raster::stack) %>%
  purrr::map(raster::unstack) %>%
  purrr::flatten() %>%
  rlang::set_names(control$bands)

# Force all layers to the same range of values
# - from OTB documentation: "Beware of potential imbalance between bands 
#   ranges as it may alter euclidean distance."
# - this also stretches spectral bands (enhances contrast)
img_stack_stretch <- img_stack %>%
  purrr::map(raster_rescale, to = c(0, 255)) %>%
  raster::stack()
raster::writeRaster(img_stack_stretch, temp_files$img_stack, overwrite = TRUE)


### Segmentation parameter pre-estimation ======================================

# Ming et al. 2015 method
# - uses spatial statistics to determine a priori best parameter values
# - largely based on average local variance (ALV)
# - designed for use with LSMS segmentation algorithm
# - for more on choosing LSMS parameters see De Luca et al. 2019

# 1. Spatial range selection

# Subset the actual site delination for analysis
img_stack_crop <- img_stack_stretch %>%
  #raster::crop(sf::st_as_sf(sf::st_buffer(tower, 75))) %>%
  raster::crop(delin) %>%
  raster::mask(delin) %>%
  raster::unstack()

# Convert sequence of spatial range values to spatial windows
ws <- seq(1, 15) * 2 + 1
# Calculate ALV at different spatial windows for each band
band_alv <- ws %>%
  purrr::map(~ purrr::map_dbl(img_stack_crop, alv, w = .x)) %>%
  purrr::map(rlang::set_names, stringr::str_c("band_", control$bands)) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(ws = ws, .before = 1)

# Assess dynamics of ALV along spatial range
band_alv <- band_alv %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    # Convert ws back to hs
    hs = floor(ws / 2),
    alv = mean(dplyr::c_across(dplyr::starts_with("band"))), .after = 1
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    # First- and second-order rate of change in ALV (Eqs. 7-8)
    # - this works best if I only evaluate based on dalv
    dalv = (alv - lag(alv)) / lag(alv),
    dalv2 = lag(dalv) - dalv
  )
# Determine spatial range parameter
hs_opt <- band_alv %>% dplyr::filter(dalv < 0.01) %>% purrr::pluck("hs", 1)
hs_opt

band_alv %>%
  dplyr::select(hs, dplyr::starts_with("band")) %>%
  tidyr::pivot_longer(-hs, names_to = "band") %>%
  ggplot2::ggplot(ggplot2::aes(hs, value, color = band)) +
  ggplot2::geom_line() +
  ggplot2::geom_point()

# 2. Spectral range selection

lv_xy <- img_stack_crop %>%
  #raster::unstack() %>%
  purrr::map(lv, w = hs_opt * 2 + 1) %>%
  raster::stack() %>%
  raster::calc(mean)
raster::plot(lv_xy)
MASS::truehist(lv_xy[])

# Optimal ranger is first peak of the histogram
lv_density <- lv_xy %>% 
  raster::values() %>% 
  na.omit() %>% 
  # Increase 'adjust' to smooth more (remove pseudo-peaks)
  density(adjust = 2)

# How many peaks?
plot(lv_density)
length(lv_density$x[which(diff(sign(diff(lv_density$y))) == -2)])

hr_opt <- lv_density %>% 
  magrittr::extract(c("x", "y")) %>%
  tibble::as_tibble() %>% 
  dplyr::mutate(d1sign = sign(y - lag(y)), d2 = d1sign - lag(d1sign)) %>% 
  # Identify histogram peaks, remove short ones
  dplyr::filter(d2 == -2, y > 0.01) %>% 
  # Select the first peak if there are multiple
  dplyr::slice_min(x) %>%
  purrr::pluck("x", 1) %>%
  sqrt() %>%
  # Parameter should be an integer in practical application
  round(0)
hr_opt

# Write segmentation parameters to file
segm_params <- append(control, list(ranger = hr_opt, spatialr = hs_opt))
sink(files_out$params)
print(segm_params) 
sink()


# Perform segmentation

# 1. Mean shift
system(paste0(
  "~/OTB-7.1.0-Darwin64/bin/", 
  "otbcli_MeanShiftSmoothing",
  # File paths
  " -in ", temp_files$img_stack,
  " -fout ", temp_files$filter_range,
  " -foutpos ", temp_files$filter_spatial,
  # Parameters
  " -ranger ", segm_params$ranger,
  " -spatialr ", segm_params$spatialr,
  " -maxiter ", segm_params$maxiter,
  " -ram ", 3072,
  " -modesearch 0"
), ignore.stdout = TRUE, ignore.stderr = TRUE)

# 2. LSMS segmentation
system(paste0(
  "~/OTB-7.1.0-Darwin64/bin/", 
  "otbcli_LSMSSegmentation",
  # File paths
  " -in ", temp_files$filter_range,
  " -inpos ", temp_files$filter_spatial,
  " -out  ", temp_files$segm_init, " uint32",
  # Parameters
  " -ranger ", segm_params$ranger,
  " -spatialr ", segm_params$spatialr,
  " -minsize 0",
  " -tilesizex ", 2500,
  " -tilesizey ", 2500
), ignore.stdout = TRUE, ignore.stderr = TRUE)

# 3. Region merging
system(paste0(
  "~/OTB-7.1.0-Darwin64/bin/",
  "otbcli_LSMSSmallRegionsMerging",
  # File paths
  " -in ", temp_files$filter_range,
  " -inseg ", temp_files$segm_init,
  " -out ", temp_files$segm_merge, " uint32",
  # Parameters
  #" -minsize ", 100,
  " -minsize ", segm_params$minsize,
  " -tilesizex ", 2500,
  " -tilesizey ", 2500, 
  " -ram ", 3072
), ignore.stdout = TRUE, ignore.stderr = TRUE)

# Get the segmented raster, check results
segments_st <- temp_files$segm_merge %>%
  stars::read_stars() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::select(segm = 1)
ggplot2::ggplot() +
  stars::geom_stars(
    data = dplyr::mutate(
      segments_st, segm = factor(stringr::str_sub(as.character(segm), -1))
    )
  ) +
  ggplot2::geom_sf(
    data = c(md$x_utm, md$y_utm) %>%
      sf::st_point() %>% 
      sf::st_sfc(crs = spatial_ref$epsg), 
    size = 2
  ) +
  ggplot2::scale_fill_brewer(palette = "Paired", guide = "none")
# Number of segments within site area
n_segm <- segments_st %>% 
  sf::st_crop(delin) %>%
  tibble::as_tibble() %>%
  tidyr::drop_na() %>%
  dplyr::pull(segm) %>%
  dplyr::n_distinct()
n_segm
as.numeric(sf::st_area(delin)) / n_segm # average size

# Write raster output
stars::write_stars(segments_st, files_out$rst)

# Also write the initial (unmerged) segments
temp_files$segm_init %>%
  stars::read_stars() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::select(segm = 1) %>%
  stars::write_stars(files_out$rst_init)

# Convert segments to polygons, write to file
segments_poly <- segments_st %>%
  sf::st_as_sf(as_points = FALSE, merge = TRUE) %>%
  dplyr::rename(segment = 1)
sf::write_sf(segments_poly, files_out$poly)



# Run initial segmentation
# segments <- SegOptim::segmentation_OTB_LSMS2(
#   # Input raster with features/bands to segment
#   inputRstPath = file.path(tempdir(), "img_stack.tif"),
#   # Algorithm parameters (OTB defaults: 15, 5, 0, 10)
#   # - over-segmentation results in greater classification accuracy (Smith 2010)
#   MS_SpectralRange = 30,
#   MS_SpatialRange = 10,
#   LSS_SpectralRange = 15,
#   LSS_SpatialRange = 15,
#   MinSize = 5,
#   lsms_maxiter = 100,
#   # Output
#   outputSegmRst = file.path(tempdir(), "segments.tif"),
#   verbose = TRUE,
#   otbBinPath = "~/OTB-7.1.0-Darwin64/bin"
# )


### Ming et al. 2012 method ====================================================

# Subset segments within site delineation
segments_eval <- segments_rst %>%
  stars::st_as_stars() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  sf::st_as_sf(as_points = FALSE, merge = TRUE) %>% 
  dplyr::rename(segment = 1) %>%
  dplyr::filter(sf::st_covered_by(geometry, delin, sparse = FALSE)) %>% 
  sf::st_join(img_bands, ., what = "inner", as_points = TRUE)

# Average band values in each segment
segments_avg <- segments_eval %>% 
  dplyr::select(-x, -y) %>% 
  dplyr::group_by(segment, geometry) %>% 
  dplyr::summarize(dplyr::across(.fns = mean), .groups = "drop")

# Homogeneity

intra_hom <- segments_eval %>%
  tibble::as_tibble() %>%
  select(segment, 3:6) %>%
  tidyr::pivot_longer(-segment, names_to = "band") %>%
  dplyr::group_by(band, segment) %>%
  dplyr::mutate(
    mean = mean(value),
    diff = (value - mean)^2
  ) %>%
  dplyr::summarize(
    area = dplyr::n(),
    homo = mean(diff), .groups = "drop_last"
  ) %>%
  dplyr::summarize(u = sqrt(sum(homo * area / sum(area))), .groups = "drop") %>%
  dplyr::summarize(u = mean(u))


# Heterogeneity
inter_het <- segments_avg %>% 
  dplyr::select(segment, 3:6, geometry) %>%
  dplyr::mutate(dplyr::across(2:5, list(mean = mean))) %>%
  sf::st_join(., ., join = sf::st_touches) %>%
  sf::st_drop_geometry() %>%
  dplyr::select(-dplyr::ends_with("mean.y")) %>%
  dplyr::rename_with(
    ~ stringr::str_replace(.x, ".x", "_value"), 
    c(-segment.x, -dplyr::contains("mean"))
  ) %>%
  dplyr::rename_with(~ stringr::str_remove(.x, ".x")) %>%
  dplyr::rename_with(~ stringr::str_replace(.x, ".y", "_adj")) %>%
  # Missing segment_adj means there are no adjacent segments
  tidyr::drop_na(segment_adj) %>%
  tidyr::pivot_longer(
    -dplyr::starts_with("segm"), names_to = c("band", "type"), names_sep = "_"
  ) %>%
  tidyr::pivot_wider(names_from = type, values_from = value) %>%
  dplyr::mutate(diff = (value - mean) * (adj - mean)) %>%
  dplyr::group_by(band, segment) %>%
  dplyr::summarize(
    mean = first(mean), 
    value = first(value),
    diff_sum = sum(diff), 
    w = dplyr::n(),
    .groups = "drop_last"
  ) %>%
  dplyr::summarize(
    v = (dplyr::n() * sum(diff_sum)) / (sum(w) * sum((value - mean)^2)),
    .groups = "drop"
  ) %>%
  dplyr::summarize(v = mean(v))
bind_cols(intra_hom,  inter_het) %>% as.data.frame()
# rgb 10.61664 0.4368428
# hsv 17.47986 0.3631229

