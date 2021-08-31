### ============================================================================
# Spatial analysis of flux area ================================================
### ============================================================================

# Purpose: 

# Input(s):

# Output(s):

settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  year = 2018, # four digit year
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

# TODO: move config to yaml file
config <- list(
  fp_model = "K15",
  features = c("ndvi", "ndwi", "z_max", "i_vp"),
  grid_width = 50,
  phi = 0.85,
  # Water level of flooded/dry threshold 
  flood_depth = 0
)

# Load the required packages
devtools::load_all("~/Projects/Flux/footprints")
devtools::load_all("~/Projects/Flux/dscalr")
library(progress)
library(lubridate)
library(tidyverse)

# Load reference files
source("~/Projects/Flux/towers/reference/site_metadata.R")


### Helper functions ===========================================================

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

cross_keys <- function(...) {
  
  dots <- rlang::list2(...)
  
  values <- dots %>%
    tidyr::crossing() %>%
    dplyr::summarize(value = prod(dplyr::c_across()))
}

cross_grids <- function(...) {
  
  dots <- rlang::list2(...)
  
  # Assume that all matrices have same dimensions
  dims <- dots %>% purrr::map(dim) %>% purrr::pluck(1)
  
  crossed <- dots %>%
    purrr::map(as.vector) %>%
    # Drop elements that only contain missing values
    purrr::discard(~ all(is.na(.x))) %>%
    purrr::lift_dl(stringr::str_c)(sep = "_")
  
  matrix(crossed, nrow = dims[1], ncol = dims[2])
}


### Initialize script settings & documentation =================================

# Load metadata file
md <- yaml::read_yaml(file.path("data", settings$site, "metadata.yml"))
tower_coords <- c(md$tower$coords$east, md$tower$coords$north)
elev_corr <- md$well$rel_elev

# Set tag for file names
tag <- make_tag(md$site_code, settings$year)

# Set the desired working directory in RStudio interface
# - assumes that the subdirectory structure is already present
path_flux <- file.path(
  "~/Projects/Flux/towers/data", settings$site, settings$year
)
path_spatial <- file.path("data", settings$site, "processing")
path_in <- file.path(path_flux, "05_footprint")
path_out <- file.path(path_spatial, "09_footprint_cover")

files <- list(
  # Daily water level
  wtd = file.path(
    "~/Projects/Flux/inundation/data/output",
    paste0(str_replace(settings$site, "JL", "F"), "_waterlevel_dd.csv")
  ),
  # Cover type raster
  class = file.path(path_spatial, "07_classification", "classes_manual"),
  # Relative elevation raster
  elev = file.path(path_spatial, "08_relative_elevation", "rel_elev.tif"),
  # Other spatial features
  feat = file.path(path_spatial, "03_feature_extraction"),
  delin = file.path(dirname(path_flux), "site_info", "delineation"),
  # Footprint data
  fp = file.path(
    path_in, "footprint", paste0("halfhourly_", config$fp_model, "_", tag)
  ),
  phi = file.path(
    path_in, "stats", 
    paste0("footprint_stats_", config$fp_model, "_", tag, ".csv")
  )
)


### Load required input data ===================================================

# Import water level data
data_wtd <- read_csv(files$wtd, show_col_types = FALSE)

# Import footprint phi data
data_phi <- read_csv(files$phi, show_col_types = FALSE)

# Import footprint grid (including AOI grid)
grid <- files$fp %>% 
  file.path("grid") %>%
  file.path(paste0(c("x", "y", "aoi"), ".txt")) %>% 
  map(read_matrix, trunc = 0) %>%
  set_names(c("x", "y", "aoi"))

# Import relative elevation
rel_elev <- raster::raster(files$elev)
flood_key <- c(dry = -1, wet = 1)

# Import classified image
class <- files$class %>%
  sf::read_sf() %>%
  transmute(
    class = if_else(class == "PH", subclass, class),
    class = if_else(is.na(class), "CA", class),
    class = factor(str_to_lower(class))
  )

# Make integer key for classes
class_key <- class %>% 
  as_tibble() %>% 
  distinct(class) %>% 
  transmute(class = levels(class), id = row_number()) %>%
  deframe()

# Read all other spatial features
feat <- files$feat %>%
  list.files(pattern = ".tif$", full.names = TRUE, recursive = TRUE) %>%
  str_subset(str_c("/", config$feat, ".", collapse = "|")) %>%
  map(raster::raster) %>%
  set_names(map(., names))


### Prepare input data =========================================================

# Convert classified polygons to raster format
class_rst <- class %>%
  stars::st_rasterize(dx = 0.6, dy = 0.6) %>%
  as("Raster")

# Map imagery onto same grid as footprint
class_grid <- class_rst %>% 
  snap_to_grid(grid, tower_coords) %>%
  with_matrix(~ as.integer(.x))

# z0_grid <- class_grid %>% 
#   replace(which(class_grid == 3), 0.5) %>% 
#   replace(which(class_grid == 2), 0.2) %>% 
#   replace(which(class_grid == 1), 0.01)

elev_grid <- snap_to_grid(rel_elev, grid, tower_coords)
feat_grid <- feat %>%
  map(~ snap_to_grid(.x, grid, tower_coords)) %>%
  prepend(list(rel_elev = elev_grid))

# Create key for combined cover types & flooding
classflood_key <- class_key %>% 
  map(~ .x * flood_key) %>% 
  imap(~ set_names(.x, str_c(.y, "_", names(.)))) %>%
  flatten() %>%
  simplify()

# Split site grid into objective quadrants
quad_grid <- with_matrix2(
  sign(grid$x), sign(grid$y), 
  ~ case_when(
    .x ==  1 & .y ==  1 ~ 1,
    .x == -1 & .y ==  1 ~ 2,
    .x == -1 & .y == -1 ~ 3,
    .x ==  1 & .y == -1 ~ 4
  )
)
quad_key <- c(q1 = 1, q2 = 2, q3 = 3, q4 = 4)
quad_cover <- quad_grid %>% 
  with_matrix2(
    class_grid, 
    ~ as.integer(fct_cross(factor(.x), factor(.y), keep_empty = TRUE))
  )
quadclass_key <- quad_cover %>%
  as.vector() %>%
  unique() %>%
  sort() %>%
  set_names(
    paste(names(quad_key), rep(names(class_key), each = 4), sep = "_")[.]
  )

### Retrieve footprints, calculate cover =======================================

# Initialize loop

# Get list of footprint files
fp_files <- files$fp %>% 
  file.path("footprints") %>%
  list.files(full.names = TRUE) %>%
  enframe(name = NULL, value = "file") %>%
  # Parse timestamps from file names
  mutate(
    timestamp = ymd_hms(str_sub(basename(file), 1, -5)),
    .before = 1
  ) %>%
  deframe()

fp_timestamps <- fp_files %>% names() %>% ymd_hms()

n_fp <- length(fp_files)

# Empty list to hold cover data
cover <- vctrs::vec_init(list(), n_fp)

# Select water level data 
wtd <- data_wtd %>%
  right_join(
    fp_timestamps %>%
      enframe(name = NULL, value = "timestamp") %>%
      mutate(date = date(timestamp)),
    by = "date"
  ) %>%
  pull(wtd_f)
p <- progress_info(n_fp)

for (i in 1:n_fp) {
  
  # Read and scale footprint matrix
  fp_temp <- read_matrix(fp_files[i], trunc = 9)

  # Mask integrated footprint to account for rapid expansion approaching 100%
  # - recommended: between 80% and 90% (Kljun et al. 2015)
  # - this means data should be filtered for phi >= p
  fp_mask <- mask_source_area(fp_temp, p = config$phi, mask_value = 0)
  fp_mask <- fp_mask * grid$aoi
  
  # Normalize integrated footprint to 1 (Tuovinen et al. 2019)
  # - the masked footprint is thus considered the 100% analytical footprint
  phi_sum <- sum(fp_mask)
  fp_norm <- fp_mask / phi_sum
  
  # Calculate water level over AOI
  # - as depth relative to land surface (negative = water below surface)
  # - give wtd a placeholder value if NA so that other covers can be summed
  aoi_wtd <- elev_corr + replace_na(wtd[i], 0.1) - elev_grid
  # Classify flooded area (1 = wet, 0 = dry)
  aoi_flooded <- with_matrix(
    aoi_wtd, ~ if_else(.x >= config$flood_depth, 1, 0)
  ) 
  
  # For cover classes, flooded is classified as (1 = wet, -1 = dry)
  aoi_cover <- class_grid * sign(aoi_wtd)
  
  # Calculate cover-type weights
  cover_wt <- summarize_cover(fp_norm, aoi_cover, levels = classflood_key)
  
  # Calculate quadrant weights
  quad_wt <- summarize_cover(fp_norm, quad_cover, levels = quadclass_key)
  
  # Calculate footprint-weighted WTD & flooded area
  wtd_wt <- summarize_cover(fp_norm, aoi_wtd, type = "numeric")
  flood_wt <- summarize_cover(fp_norm, aoi_flooded, type = "numeric")
  feat_wt <- feat_grid %>%
    map(~ summarize_cover(fp_norm, .x, type = "numeric")) %>%
    simplify()
  
  # Add all weights to list 
  cover[[i]] <- c(
    phi_mask = phi_sum, cover_wt, "wtd_wt" = wtd_wt, "flood_wt" = flood_wt, 
    feat_wt, quad_wt
  )

  p$tick()
}

# Gather everything into one data frame
data_cover <- cover %>%
  bind_rows() %>%
  mutate(timestamp = fp_timestamps) %>%
  relocate(timestamp) %>%
  right_join(select(data_phi, timestamp), by = "timestamp") %>%
  arrange(timestamp)

# Sum crosstabs
data_cover_patch <- data_cover %>% 
  select(timestamp, ends_with(c("wet", "dry"))) %>% 
  pivot_longer(-timestamp, names_to = c("patch", "flood"), names_sep = "_") %>% 
  group_by(timestamp, patch) %>% 
  summarize(p = sum(value), .groups = "drop") %>% 
  pivot_wider(names_from = patch, values_from = p) %>%
  rowwise() %>%
  mutate(
    # 3-patch scheme (open water, herbaceous vegetation, forested wetland)
    wtr = sum(c_across(any_of(c("av", "ow", "pe")))),
    veg = sum(c_across(any_of(c("ca", "el", "ju", "pa", "hv")))),
    frs = sum(c_across(any_of(c("fo", "fw")))),
    # Sedges/rushes
    sr = sum(c_across(any_of(c("ca", "el", "ju")))),
    .after = 2
  ) %>%
  ungroup() %>%
  bind_cols(select(data_cover, phi_mask))

data_cover_flood <- data_cover %>%
  mutate(date = date(timestamp)) %>%
  left_join(select(data_wtd, date, wtd = wtd_f), by = "date") %>%
  select(timestamp, phi_mask, wtd, ends_with(c("wet", "dry"))) %>%
  # Set wet/dry covers to NA if WTD was missing
  mutate(
    across(-c(timestamp, phi_mask), ~ if_else(is.na(wtd), NA_real_, .x))
  ) %>%
  rowwise() %>%
  mutate(
    wet = sum(c_across(ends_with("wet"))),
    dry = sum(c_across(ends_with("dry")))
  ) %>%
  mutate(
    # 3-patch scheme (open water, herbaceous vegetation, forested wetland)
    wtr_wet = sum(c_across(any_of(paste0(c("av", "ow", "pe"), "_wet")))),
    wtr_dry = sum(c_across(any_of(paste0(c("av", "ow", "pe"), "_dry")))),
    veg_wet = sum(c_across(any_of(paste0(c("ca", "el", "ju", "pa"), "_wet")))),
    veg_dry = sum(c_across(any_of(paste0(c("ca", "el", "ju", "pa"), "_dry")))),
    frs_wet = sum(c_across(any_of("fo_wet"))),
    frs_dry = sum(c_across(any_of("fo_dry"))),
    # Sedges/rushes
    sr_wet = sum(c_across(any_of(paste0(c("ca", "el", "ju"), "_wet")))),
    sr_dry = sum(c_across(any_of(paste0(c("ca", "el", "ju"), "_dry"))))
  ) %>%
  ungroup()
  
data_cover_quad <- data_cover %>%
  select(timestamp, phi_mask, starts_with(c("q1", "q2", "q3", "q4"))) %>%
  rowwise() %>%
  mutate(
    q1 = sum(c_across(starts_with("q1"))),
    q2 = sum(c_across(starts_with("q2"))),
    q3 = sum(c_across(starts_with("q3"))),
    q4 = sum(c_across(starts_with("q4")))
  ) %>%
  mutate(
    # 3-patch scheme (open water, herbaceous vegetation, forested wetland)
    q1_wtr = sum(c_across(any_of(paste0("q1_", c("av", "ow", "pe"))))),
    q2_wtr = sum(c_across(any_of(paste0("q2_", c("av", "ow", "pe"))))),
    q3_wtr = sum(c_across(any_of(paste0("q3_", c("av", "ow", "pe"))))),
    q4_wtr = sum(c_across(any_of(paste0("q4_", c("av", "ow", "pe"))))),
    q1_veg = sum(c_across(any_of(paste0("q1_", c("ca", "el", "ju", "pa"))))),
    q2_veg = sum(c_across(any_of(paste0("q2_", c("ca", "el", "ju", "pa"))))),
    q3_veg = sum(c_across(any_of(paste0("q3_", c("ca", "el", "ju", "pa"))))),
    q4_veg = sum(c_across(any_of(paste0("q4_", c("ca", "el", "ju", "pa"))))),
    q1_frs = sum(c_across(any_of("q1_fo"))),
    q2_frs = sum(c_across(any_of("q2_fo"))),
    q3_frs = sum(c_across(any_of("q3_fo"))),
    q4_frs = sum(c_across(any_of("q4_fo"))),
    # Sedges/rushes
    q1_sr = sum(c_across(any_of(paste0("q1_", c("ca", "el", "ju"))))),
    q2_sr = sum(c_across(any_of(paste0("q2_", c("ca", "el", "ju"))))),
    q3_sr = sum(c_across(any_of(paste0("q3_", c("ca", "el", "ju"))))),
    q4_sr = sum(c_across(any_of(paste0("q4_", c("ca", "el", "ju")))))
  ) %>%
  ungroup()

# Peek distributions
# data_cover_all %>%
#   pivot_longer(c(wet, dry)) %>%
#   ggplot(aes(value, fill = name)) +
#   geom_density(na.rm = TRUE, alpha = 0.6)
# 
# data_cover_all %>%
#   summarize(across(c(-timestamp, -phi_mask), sum, na.rm = TRUE)) %>%
#   pivot_longer(everything())

# Quadrant weights
# ggplot() + 
#   stars::geom_stars(
#     data = quad_grid %>% 
#       apply(2, rev) %>% t() %>%
#       stars::st_as_stars() %>% 
#       as_tibble(center = FALSE) %>% 
#       left_join(
#         data_cover_all %>% 
#           filter(phi_mask > 0.8) %>%
#           select(q1, q2, q3, q4) %>% 
#           summarize(across(.fns = ~ sum(.x, na.rm = TRUE))) %>% 
#           pivot_longer(everything()) %>%
#           mutate(name = as.numeric(str_remove_all(name, "q"))),
#         by = c("A1" = "name")
#       ) %>% 
#       select(-A1) %>% 
#       stars::st_as_stars()
#   ) + 
#   scale_fill_distiller(palette = "Spectral", trans = "log") + 
#   coord_fixed()

# Write output to file
patch_out <- file.path(path_out, paste0("footprint_cover_patch_", tag, ".csv"))
flood_out <- file.path(path_out, paste0("footprint_cover_flood_", tag, ".csv"))
quad_out <- file.path(path_out, paste0("footprint_cover_quad_", tag, ".csv"))

write_csv(data_cover_patch, patch_out)
write_csv(data_cover_flood, flood_out)
write_csv(data_cover_quad, quad_out)


