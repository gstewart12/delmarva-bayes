
rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

# References

# Li, J., Tran, M., & Siwabessy, J. (2016). Selecting Optimal Random Forest 
# Predictive Models: A Case Study on Predicting the Spatial Distribution of 
# Seabed Hardness. PLOS ONE, 11(2), e0149089. 
# https://doi.org/10.1371/journal.pone.0149089
# 
# Ma, L., Fu, T., Blaschke, T., Li, M., Tiede, D., Zhou, Z., et al. (2017). 
# Evaluation of Feature Selection Methods for Object-Based Land Cover Mapping of 
# Unmanned Aerial Vehicle Imagery Using Random Forest and Support Vector Machine 
# Classifiers. ISPRS International Journal of Geo-Information, 6(2), 51. 
# https://doi.org/10.3390/ijgi6020051
# 
# Mikola, J., Virtanen, T., Linkosalmi, M., Vähä, E., Nyman, J., Postanogova, 
# O., et al. (2018). Spatial variation and linkages of soil and vegetation in 
# the Siberian Arctic tundra – coupling field observations with remote sensing 
# data. Biogeosciences, 15(9), 2781–2801. 
# https://doi.org/10.5194/bg-15-2781-2018
#
# Millard, K., & Richardson, M. (2015). On the Importance of Training Data 
# Sample Selection in Random Forest Image Classification: A Case Study in 
# Peatland Ecosystem Mapping. Remote Sensing, 7(7), 8489–8515. 
# https://doi.org/10.3390/rs70708489
# 
# Räsänen, A., Kuitunen, M., Tomppo, E., & Lensu, A. (2014). Coupling 
# high-resolution satellite imagery with ALS-based canopy height model and 
# digital elevation model in object-based boreal forest habitat type 
# classification. ISPRS Journal of Photogrammetry and Remote Sensing, 94, 
# 169–182. https://doi.org/10.1016/j.isprsjprs.2014.05.003




control <- list(
  # Maximum iterations for determining variable importance (Mikola et al. 2018) 
  max_boruta_runs = 1000,
  # Training split proportion (Millard & Richardson 2015)
  split_prop = 0.795,
  # Cutoff for removing highly correlated predictors (Li et al. 2016)
  corr_thr = 0.99,
  # Number of RF trees to grow (Millard & Richardson 2015)
  n_trees = 1000,
  n_ens_models = 100
)

#library(SegOptim)
library(progress)
library(tidyverse)

raster_rescale <- function(x, ...) {
  values <- raster::values(x)
  new_values <- scales::rescale(values, ...)
  raster::setValues(x, new_values)
}

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
  seeds = "/Users/Graham/Desktop/DATA/Flux/tools/reference/seeds.R",
  spatial_ref = "/Users/Graham/Desktop/DATA/Flux/tools/reference/spatial.R",
  # Path to site area
  delin = file.path(dirname(dirname(wd)), "site_info", "delineation"),
  point = "/Users/Graham/Desktop/DATA/Spatial/survey/flux_pts",
  # Path to image segments
  segm = file.path(wd, "04_segmentation", "segments_rgbn.tif"),
  # Path to segment features
  feat = file.path(wd, "05_segment_features", "segment_features.csv"),
  # Path to training points
  train = file.path(wd, "06_training_points", "training_points"),
  # Path to output files
  out = file.path(wd, "07_classification")
)

# Load reference files
source(paths$seeds)
source(paths$spatial_ref)
seeds <- purrr::pluck(seeds, settings$site, "obia")


### Load input data ============================================================

# Read site area polygon
delin <- sf::read_sf(paths$delin)

# Read points of interest
point <- sf::read_sf(paths$point)
tower <- point %>% 
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::filter(site == settings$site, type == "tower") %>%
  sf::st_geometry()

# Load image segments
segm <- paths$segm %>% 
  stars::read_stars() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::select(segment = 1)

# Get train point data
point_train <- sf::read_sf(paths$train)

# Prep training data
point_train <- point_train %>%
  dplyr::select(-dplyr::any_of("segment"), -certainty, -subclass) %>%
  dplyr::mutate(
    ID = as.integer(ID),
    class = as.integer(factor(class))
  ) %>%
  sf::st_join(sf::st_as_sf(segm, as_points = FALSE)) %>%
  # Keep only one point per segment
  #dplyr::distinct(segment, .keep_all = TRUE) %>%
  dplyr::select(-segment, -ID)

# Import segment feature data
data_feat <- readr::read_csv(
  paths$feat, col_types = readr::cols(.default = readr::col_guess()), 
  progress = FALSE
)

# Convert raster data to tidy format
# - this is a hacky way of doing it but c.stars simply doesn't work (why??) 
# feat_st <- feat_rst %>% 
#   stars::st_as_stars() %>% 
#   tibble::as_tibble(center = FALSE) %>% 
#   dplyr::rename(value = 4) %>% 
#   tidyr::pivot_wider(names_from = band, values_from = value) %>% 
#   dplyr::rename_with(~ stringr::str_replace(.x, "\\.", "_")) %>%
#   stars::st_as_stars() %>%
#   sf::st_set_crs(spatial_ref$epsg)


### Assign features to segments ================================================

# Convert segments to polygon format
segm_poly <- sf::st_as_sf(segm, as_points = FALSE, merge = TRUE)

# Assign training classes to segments
segm_class <- segm_poly %>%
  sf::st_join(., point_train, what = "inner", as_points = FALSE) %>%
  dplyr::filter(!is.na(class)) %>% 
  # Remove duplicated segments
  # - happens if multiple training points end up within one segment
  dplyr::distinct(segment, .keep_all = TRUE) %>%
  dplyr::mutate(class = factor(class))

# Check training segments
segm_class %>%
  ggplot2::ggplot() + 
  ggplot2::geom_sf(ggplot2::aes(fill = class, color = class))

# Class balance
table(segm_class$class)

# Add point classifications to feature data
data_feat <- dplyr::left_join(
  data_feat, sf::st_drop_geometry(segm_class), by = "segment"
)

# Subset classified training segments
data_class <- data_feat %>%
  dplyr::relocate(class) %>%
  dplyr::filter(!is.na(class)) %>%
  dplyr::mutate(class = factor(class))


### Feature selection ==========================================================

# Necessity & implementation described in Ma et al. 2017
# - important especially for OBIA due to increased feature space

# 1. Remove features that are highly correlated to other features
# - done before variable importance (Li et al. 2016)
keep_corr <- data_class %>% 
  dplyr::select(-segment, -class) %>%
  recipes::recipe(~ .) %>% 
  # Spearman rank method since variables generally not normally distributed
  recipes::step_corr(
    dplyr::everything(), threshold = control$corr_thr, method = "spearman"
  ) %>%
  recipes::prep() %>%
  recipes::juice() %>%
  names()
keep_corr

# Select relevant features to create training set
data_class_corr <- dplyr::select(
  data_class, class, segment, dplyr::all_of(keep_corr)
)

# 2. Remove features deemed unimportant for classification
# - tends to improve classification performance (Rasanen et al. 2014)
set.seed(seeds$boruta)
# - ~5 mins for ~200x600 df
var_imp <- Boruta::Boruta(
  class ~ ., 
  data = dplyr::select(data_class_corr, -segment), 
  maxRuns = control$max_boruta_runs
)
(var_imp <- Boruta::TentativeRoughFix(var_imp))
(keep_imp <- Boruta::getSelectedAttributes(var_imp))

# Variable importance
keep_var_imp <- var_imp %>% 
  purrr::pluck("ImpHistory") %>% 
  tibble::as_tibble() %>% 
  dplyr::select(-dplyr::contains("Shadow")) %>%
  dplyr::summarize(dplyr::across(.fns = mean)) %>%
  tidyr::pivot_longer(dplyr::everything()) %>% 
  dplyr::filter(is.finite(value)) %>%
  dplyr::arrange(dplyr::desc(value))
keep_var_imp

# Select relevant features to create training set
data_class_imp <- dplyr::select(
  data_class_corr, class, segment, dplyr::all_of(keep_imp)
)

# Write to file
readr::write_csv(data_class_imp, file.path(paths$out, "selected_features.csv"))
# Read back in (so feature selection can be skipped)
data_class_imp <- readr::read_csv(
  file.path(paths$out, "selected_features.csv"), 
  col_types = readr::cols(
    class = readr::col_factor(levels = c(1, 2, 3)), 
    .default = readr::col_double()
  )
)

# Which vars remain?
keep_imp %>%
  tibble::enframe(name = NULL, value = "name") %>%
  tidyr::separate(
    name, c("var", "desc", "stat"), extra = "merge", fill = "left"
  ) %>% 
  tidyr::unite("var", 1:2, na.rm = TRUE) %>% 
  dplyr::count(var, sort = TRUE)


### Classification =============================================================

# Set preprocessing recipe
rf_rec <- data_class_imp %>% 
  recipes::recipe(class ~ .) %>% 
  recipes::update_role(segment, new_role = "ID") %>%
  recipes::step_range(recipes::all_predictors(), min = 0, max = 1)

# Set model
rf_mod <- parsnip::rand_forest() %>%
  parsnip::set_args(trees = control$n_trees) %>%
  parsnip::set_mode("classification") %>%
  parsnip::set_engine("randomForest")

# Set model fit workflow
rf_wflow <- workflows::workflow() %>%
  workflows::add_model(rf_mod) %>% 
  workflows::add_recipe(rf_rec) 

# First run: using all data for training

# Fit model
set.seed(seeds$fit)
fit_all <- parsnip::fit(rf_wflow, data = data_class_imp)

# Get predictions
pred_segm_all <- segm_poly %>%
  dplyr::right_join(dplyr::select(data_feat, segment), ., by = "segment") %>%
  dplyr::bind_cols(
    predict(fit_all, new_data = data_feat),
    predict(fit_all, new_data = data_feat, type = "prob")
  ) %>%
  sf::st_as_sf()

pred_all <- pred_segm_all %>%
  sf::st_drop_geometry() %>%
  tibble::as_tibble() %>%
  dplyr::right_join(tibble::as_tibble(segm, center = FALSE), by = "segment") %>%
  dplyr::relocate(x, y) %>%
  stars::st_as_stars() %>%
  sf::st_set_crs(spatial_ref$epsg)

# Check predictions
pred_all %>%
  sf::st_crop(sf::st_bbox(delin), as_points = FALSE) %>%
  plot_stars(.pred_class) +
  ggplot2::geom_sf(data = delin, fill = NA) +
  ggplot2::theme(legend.position = "none")

# Write predicted classes as raster file
class_pred_all <- dplyr::select(pred_all, class = .pred_class)
stars::write_stars(class_pred_all, file.path(paths$out, "classes_all.tif"))
class_prob_all <- pred_all %>%
  dplyr::select(prob_1 = .pred_1, prob_2 = .pred_2, prob_3 = .pred_3) %>%
  stars::st_redimension(along = list(band = names(.)))
stars::write_stars(class_prob_all, file.path(paths$out, "classes_prob_all.tif"))


# Second run: single model

# Split into training (80%) and testing sets (20%)
# - for some reason need to set prop to 0.795 to get 160/40 split
set.seed(seeds$split)
data_split <- rsample::initial_split(data_class_imp, prop = control$split_prop)
data_train <- rsample::training(data_split)
table(data_train$class)
data_test  <- rsample::testing(data_split)
table(data_test$class)

# Fit model
set.seed(seeds$fit)
fit_sgl <- parsnip::fit(rf_wflow, data = data_train)

# Model performance - independent evaluation
fit_sgl %>% 
  predict(data_test) %>%
  dplyr::bind_cols(dplyr::select(data_test, class)) %>%
  yardstick::accuracy(class, .pred_class)
# Model performance - overall
fit_sgl %>% 
  predict(data_class_imp) %>%
  dplyr::bind_cols(dplyr::select(data_class_imp, class)) %>%
  yardstick::accuracy(class, .pred_class)

# Variable importance
fit_sgl %>% 
  workflows::pull_workflow_fit() %>% 
  vip::vip(num_features = 20)

# Get predictions
pred_segm_sgl <- segm_poly %>%
  dplyr::right_join(dplyr::select(data_feat, segment), ., by = "segment") %>%
  dplyr::bind_cols(
    predict(fit_sgl, new_data = data_feat),
    predict(fit_sgl, new_data = data_feat, type = "prob")
  ) %>%
  sf::st_as_sf()

pred_sgl <- pred_segm_sgl %>%
  sf::st_drop_geometry() %>%
  tibble::as_tibble() %>%
  dplyr::right_join(tibble::as_tibble(segm, center = FALSE), by = "segment") %>%
  dplyr::relocate(x, y) %>%
  stars::st_as_stars() %>%
  sf::st_set_crs(spatial_ref$epsg)

# Check predictions
pred_sgl %>%
  sf::st_crop(sf::st_bbox(delin), as_points = FALSE) %>%
  plot_stars(.pred_class) +
  ggplot2::geom_sf(data = delin, fill = NA) +
  ggplot2::theme(legend.position = "none")

# Write predicted classes as raster file
class_pred_sgl <- dplyr::select(pred_sgl, class = .pred_class)
stars::write_stars(class_pred_sgl, file.path(paths$out, "classes_sgl.tif"))
class_prob_sgl <- pred_sgl %>%
  dplyr::select(prob_1 = .pred_1, prob_2 = .pred_2, prob_3 = .pred_3) %>%
  stars::st_redimension(along = list(band = names(.)))
stars::write_stars(class_prob_sgl, file.path(paths$out, "classes_prob_sgl.tif"))


# Second run: model ensemble

# Resample to desired number of train/test splits
set.seed(seeds$split)
data_res <- seq(1, control$n_ens_models) %>%
  purrr::map(
    ~ rsample::validation_split(data_class_imp, prop = control$split_prop)
  ) %>% 
  dplyr::bind_rows() %>%
  dplyr::mutate(id = stringr::str_c(id, dplyr::row_number()))

# Fit the models
set.seed(seeds$fit)
fit_res <- data_res %>%
  dplyr::mutate(
    train = purrr::map(splits, rsample::analysis),
    test = purrr::map(splits, rsample::assessment),
    fit = purrr::map(train, ~ parsnip::fit(rf_wflow, data = .x)),
    trees = rf_wflow %>%
      workflows::pull_workflow_spec() %>% 
      purrr::pluck("args", "trees") %>% 
      rlang::eval_tidy(),
    oob = fit %>% 
      purrr::map(tune::extract_model) %>%
      purrr::map(purrr::pluck, "err.rate") %>%
      purrr::map2_dbl(trees, ~ purrr::pluck(.x, .y)) %>%
      magrittr::subtract(1, .),
    pred = fit %>% 
      purrr::map2(test, ~ predict(.x, .y)) %>%
      purrr::map2(test, ~ dplyr::bind_cols(.x, dplyr::select(.y, class))),
    accuracy = pred %>%
      purrr::map(yardstick::accuracy, class, .pred_class) %>%
      purrr::map_dbl(dplyr::pull, .estimate)
  )

# Model performance - independent evaluation
dplyr::summarize(fit_res, dplyr::across(c(oob, accuracy), mean))

# Variable importance
# TODO write to a file somewhere
vi_res <- fit_res %>%
  dplyr::transmute(
    model = stringr::str_replace_all(id, "validation", "v"),
    vi = fit %>% 
      purrr::map(workflows::pull_workflow_fit) %>%
      purrr::map(vip::vi) %>%
      purrr::map(dplyr::rename, var = 1, imp = 2) %>%
      purrr::map(tibble::rowid_to_column, var = "rank")
  ) 
vi_res %>%
  tidyr::unnest(vi) %>%
  dplyr::group_by(var) %>%
  dplyr::summarize(dplyr::across(c(rank, imp), mean), .groups = "drop") %>%
  dplyr::arrange(rank)

# Get model predictions
pred_res <- fit_res %>%
  #dplyr::arrange(dplyr::desc(accuracy), dplyr::desc(oob)) %>%
  #dplyr::slice_head(n = 50) %>%
  dplyr::transmute(
    model = stringr::str_replace_all(id, "validation", "v"),
    segment = data_feat %>% 
      dplyr::select(segment) %>% 
      rlang::list2(),
    pred_c = purrr::map(fit, predict, new_data = data_feat),
    pred_p = purrr::map(fit, predict, new_data = data_feat, type = "prob"),
    pred = purrr::pmap(list(segment, pred_c, pred_p), dplyr::bind_cols)
  ) %>%
  dplyr::select(model, pred)

# Get ensemble predictions
pred_ens_segm <- pred_res %>% 
  tidyr::unnest(pred) %>% 
  dplyr::group_by(segment) %>% 
  dplyr::summarize(
    count = list(vctrs::vec_count(.pred_class)), 
    class = count %>% 
      purrr::map_int(purrr::pluck, "key", 1) %>% 
      forcats::as_factor(),
    n = dplyr::n(),
    dplyr::across(c(.pred_3, .pred_2, .pred_1), mean),
    .groups = "drop"
  ) %>% 
  tidyr::unnest(count) %>%
  tidyr::pivot_wider(
    names_from = key, names_glue = "prob_{key}", values_from = count
  ) %>%
  dplyr::mutate(
    class = forcats::fct_inseq(class),
    dplyr::across(prob_3:prob_1, ~ tidyr::replace_na(.x / n, 0))
  ) %>%
  dplyr::select(-n)

pred_ens <- pred_ens_segm %>%
  dplyr::right_join(tibble::as_tibble(segm, center = FALSE), by = "segment") %>%
  dplyr::relocate(x, y) %>%
  stars::st_as_stars() %>%
  sf::st_set_crs(spatial_ref$epsg)

# Model performance - overall
pred_ens_segm %>% 
  dplyr::select(segment, pred_class = class) %>%
  dplyr::right_join(
    dplyr::select(data_class_imp, segment, class), by = "segment"
  ) %>%
  yardstick::accuracy(class, pred_class)

# Check ensemble predictions
pred_ens %>%
  #sf::st_crop(sf::st_bbox(delin), as_points = FALSE) %>%
  #dplyr::mutate(class = factor(class)) %>%
  plot_stars(class) +
  ggplot2::geom_sf(data = delin, fill = NA) +
  ggplot2::theme(legend.position = "none")

# Check ensemble stability
pred_ens %>%
  #sf::st_crop(sf::st_bbox(delin), as_points = FALSE) %>%
  dplyr::select(prob_3:prob_1) %>%
  stars::st_redimension(along = list(prob = names(.))) %>%
  plot_stars() +
  ggplot2::facet_wrap(~ prob) +
  ggplot2::geom_sf(data = delin, fill = NA) +
  ggplot2::scale_fill_distiller(
    palette = "GnBu", direction = 1, trans = "log1p"
  ) +
  ggplot2::theme(legend.position = "none")

# Write predicted classes as raster file
class_pred_ens <- dplyr::select(pred_ens, class)
stars::write_stars(class_pred_ens, file.path(paths$out, "classes_ens.tif"))
class_prob_ens <- pred_ens %>%
  dplyr::select(prob_1, prob_2, prob_3) %>%
  stars::st_redimension(along = list(band = names(.)))
stars::write_stars(class_prob_ens, file.path(paths$out, "classes_prob_ens.tif"))

# Write predicted classes as shapefile
class_pred_segm <- segm_poly %>% 
  dplyr::left_join(pred_ens_segm, by = "segment") %>% 
  #dplyr::select(-dplyr::starts_with(".pred")) %>%
  dplyr::mutate(class_orig = class, .after = 5)
if (!dir.exists(file.path(paths$out, "classes_segm"))) {
  dir.create(file.path(paths$out, "classes_segm"))
}
sf::write_sf(
  class_pred_segm, file.path(paths$out, "classes_segm", "classes_segm.shp")
)
# Make a copy for revision, if necessary
if (!dir.exists(file.path(paths$out, "classes_segm_rev"))) {
  dir.create(file.path(paths$out, "classes_segm_rev"))
  sf::write_sf(
    class_pred_segm, 
    file.path(paths$out, "classes_segm_rev", "classes_segm_rev.shp")
  )
}

# (Load the revised shapefile back here to convert to raster)
class_pred_segm_rev <- sf::read_sf(file.path(paths$out, "classes_segm_rev"))
class_pred_rev <- class_pred_segm_rev %>%
  sf::st_drop_geometry() %>%
  tibble::as_tibble() %>%
  dplyr::right_join(tibble::as_tibble(segm, center = FALSE), by = "segment") %>%
  dplyr::mutate(dplyr::across(dplyr::contains("class"), as.integer)) %>%
  dplyr::select(x, y, class) %>%
  stars::st_as_stars() %>%
  sf::st_set_crs(spatial_ref$epsg)
stars::write_stars(class_pred_rev, file.path(paths$out, "classes_rev.tif"))



###
# Set new patches within classified image

spectral_distance <- function(x, y, p = 2, type = c("minkowski", "sa")) {
  
  # x and y are vectors of same length
  type <- rlang::arg_match(type)
  
  if (type == "minkowski") {
    # manhattan distance: p = 1; euclidean_distance: p = 2
    dist <- sqrt(sum(abs(x - y)^p))
  } else if (type == "sa") {
    # spectral angle
    dist <- acos(sum(x * y) / sqrt(sum(x^2) * sum(y^2)))
  }
  
  dist
}
min_size <- 50

class_pred_ens <- file.path(paths$out, "classes_ens.tif") %>% 
  stars::read_stars() %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  dplyr::select(class = 1)

class_pred_segm <- segm_poly %>%
  # Using segmentation features - maybe better to use more meaningful features?
  # - e.g. canopy height, DEM, ndvi
  dplyr::left_join(dplyr::select(
    data_feat, segment, hsv_1_mean, hsv_2_mean, hsv_3_mean, nir_mean
  ), by = "segment") %>%
  sf::st_join(
    sf::st_as_sf(class_pred_ens, as_points = FALSE, merge = TRUE), 
    join = sf::st_covered_by
  )

# Rescale features for accurate Euclidean distances
patches <- dplyr::mutate(
  class_pred_segm, 
  dplyr::across(c(-segment, -class, -geometry), scales::rescale)
)
p <- progress_info(nrow(patches)) # takes ~10 min to run
i <- 1 # iterator start
repeat ({
  
  if (i > nrow(patches)) {
    break
  }
  
  patch <- patches %>% 
    dplyr::slice(i) %>% 
    tibble::as_tibble() %>% 
    tidyr::nest(values = c(-segment, -geometry, -class)) %>% 
    as.list() %>% 
    purrr::map_at("values", ~ purrr::as_vector(purrr::flatten(.x)))
  
  # Large segments don't need to be merged
  # - but they are still available for joining to smaller ones
  if (as.numeric(sf::st_area(patch$geometry)) >= min_size) {
    i <- i + 1
    p$tick()
    next
  }
  
  adj <- patch$geometry %>% 
    # Shared corners & the patch itself don't count (Rook's case)
    sf::st_relate(patches, pattern = "F***1****") %>%
    purrr::pluck(1) %>% 
    dplyr::slice(patches, .)
  
  # "Island" segments cannot be merged with anything
  if (nrow(adj) == 0) {
    i <- i + 1
    p$tick()
    next
  }
  
  adj <- adj %>%
    dplyr::filter(class == patch$class) %>%
    sf::st_drop_geometry() %>% 
    dplyr::select(-class)
  
  # Find adjacent segment with lowest spectral distance
  adj_segment <- adj %>% 
    dplyr::rowwise() %>% 
    dplyr::mutate(
      dist = spectral_distance(patch$values, dplyr::c_across(-segment))
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::arrange(dist) %>% 
    purrr::pluck("segment", 1)
  
  # Update geometry
  # - this is more complicated than group_by/summarize w/ all data, but faster
  # - still very slow; there must be a better way to do this
  # - could calculate all areas beforehand, keep track of new areas by addition
  new_patch <- patches %>% 
    purrr::assign_in(
      list("segment", which(patches$segment == patch$segment)), adj_segment
    ) %>% 
    dplyr::filter(segment == adj_segment) %>% 
    dplyr::group_by(segment) %>% 
    dplyr::summarize(dplyr::across(-geometry, .fns = mean), .groups = "drop")
  patches <- patches %>%
    dplyr::filter(!segment %in% c(patch$segment, adj_segment)) %>%
    dplyr::bind_rows(new_patch)
    #dplyr::arrange(segment)
  
  p$tick()
})

# Check patches
patches %>%
  sf::st_crop(sf::st_buffer(tower, 50)) %>%
  ggplot() +
  geom_sf(aes(fill = class))

# Number of patches within 50-m tower radius
patches %>%
  sf::st_crop(sf::st_buffer(tower, 50)) %>%
  nrow()

# Avg. patch area
patches %>%
  dplyr::mutate(area = sf::st_area(geometry)) %>%
  dplyr::summarize(mean(area))

# Write patches to file
patches %>%
  dplyr::select(segment) %>%
  stars::st_rasterize(template = segm) %>%
  sf::st_set_crs(spatial_ref$epsg) %>%
  stars::write_stars(file.path(paths$out, "patches_50.tif"))

