
# Patch type flux modeling =====================================================


# Resources
# https://cchecastaldo.github.io/BayesianShortCourse/Syllabus.html

# Mixture model - fundamental concepts
# fch4 = p_wtr * fch4_wtr + p_veg * fch4_veg + p_frs * fch4_frs
# F = phi_a * f_a + phi_b * f_b + phi_c * f_c
# - have phi, need f
# - assume f is normally distributed, but cannot assume phi is
# - cannot solve algebraically, but can solve probabilistically 

# JAGS format
# model {
#   # priors
#   beta1 ~ dnorm(mu.beta1, tau.beta1)
#   beta2 ~ dnorm(mu.beta2, tau.beta2)
#   sigma ~ dunif(mu.sigma, tau.sigma)
#   tau <- 1 / sigma^2
#   
#   # likelihood
#   for (i in 1:length(y)) {
#     # mu is the "true" flux
#     mu[i] <- beta1 * x1[i] + beta2 * x2[i]
#     # y is the observed flux which is considered to be drawn from a distribution 
#   described by mu and tau
#     y[i] ~ dnorm(mu[i], tau)
#   }
# }

rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLR", # three letter site code
  years = c(2018, 2019), # vector of four digit years
  dir = "/Users/Graham/Desktop/DATA", # where all data can be found
  date = lubridate::today(), # date script was run
  info = devtools::session_info() # R session info: R, OS, packages
)

control <- list(
  flux = rlang::expr(fch4),
  # Parameters for prior (uniform) distribution of sigma
  sigma_prior_min = 0,
  sigma_prior_max = 100,
  # Number of independent chains of MCMC sampling
  n_chains = 3,
  # Early of samples to discard at the beginning of chains
  n_adapt = 1000,
  # Number of samples to take from the posterior distribution
  n_update = 10000,
  n_iter = 10000,
  pred_fun = quantile
)

devtools::load_all("/Users/Graham/R Projects/footprints")
library(tidyverse)

set_jags_inits <- function(priors, n_chains) {
  
  n_chains %>%
    vctrs::vec_init(list(), .) %>%
    purrr::map(
      ~ priors %>% 
        magrittr::extract(c("mu", "sigma")) %>% 
        purrr::pmap(function(mu, sigma) rnorm(1, mu, sigma)) %>%
        rlang::set_names(priors$param) %>%
        purrr::prepend(list(
          sigma = runif(1, priors$sigma_y[1], priors$sigma_y[2])
        ))
    )
}

set_jags_model <- function(priors, posterior_check = FALSE) {
  
  beta_priors <- priors %>% 
    magrittr::extract(c("param", "mu", "tau")) %>% 
    purrr::pmap(
      function(param, mu, tau) glue::glue("{param} ~ dnorm({mu}, {tau})")
    ) %>%
    purrr::map(stringr::str_c, "\n") %>%
    stringr::str_flatten(collapse = " ")
  
  likelihood <- priors$param %>%
    rlang::set_names() %>%
    purrr::map_chr(~ stringr::str_replace(.x, "beta", "x")) %>%
    purrr::imap_chr(~ stringr::str_c(.y, " * ", .x, "[i]")) %>%
    stringr::str_flatten(collapse = " + ")
  
  if (posterior_check) {
    # Simulated data for posterior predictive checks
    predictions <- glue::glue("
      y.sim[i] ~ dnorm(mu[i], tau) 
      sq.data[i] <- (y[i] - mu[i])^2
      sq.sim[i] <- (y.sim[i] - mu[i])^2
    ")
    
    # Bayesian p values
    # - probability that the test statistic calculated from simulated data is
    #   more extreme than that from the observed data
    p_values <- glue::glue("
      sd.data <- sd(y[])
      sd.sim <- sd(y.sim[])
      p.sd <- step(sd.sim - sd.data)
    
      mean.data <- mean(y[])
      mean.sim  <- mean(y.sim[])
      p.mean <- step(mean.sim - mean.data)
    
      discrep.data <- sum(sq.data[])
      discrep.sim <- sum(sq.sim[])
      p.discrep <- step(discrep.sim - discrep.data)
    ")
    
  } else {
    predictions <- ""
    p_values <- ""
  }
  
  model <- glue::glue("
    model {{

    # priors
    {beta_priors}
    sigma ~ dunif({priors$sigma_y[1]}, {priors$sigma_y[2]}) 
    tau <- 1 / sigma^2
  
    # likelihood
    for (i in 1:length(y)) {{
      mu[i] <- {likelihood}
      y[i] ~ dnorm(mu[i], tau)
      
      {predictions}
    }}
    
    {p_values}
  
    }}
  ")
  
  as.character(model)
}

subset.mcmc <- function(x, vars) {
  
  if (all(vars == "all")) {
    return(x)
  }
  
  x[, colnames(x[[1]]) %in% vars]
}

diagnose.mcmc <- function(x, ...) {
  
  # Gelman and Rubin diagnostic
  # - point estimates and 97.5% quantiles should approach 1
  # - if 95% quantile > 1.05, more iterations should be run 
  gelman <- coda::gelman.diag(x)
  gelman_tidy <- gelman %>%
    purrr::pluck("psrf") %>%
    tibble::as_tibble(rownames = "var") %>%
    dplyr::rename(psrf_point = 2, psrf_upper = 3) %>%
    dplyr::mutate(mpsrf = purrr::pluck(gelman, "mpsrf")) %>%
    dplyr::rename_with(~ stringr::str_c(.x, "_gelman"), -var)
  
  # Heidelberger and Welch diagnostic
  # - all chains and parameters should pass stationary & half-width tests (= 1)
  # - if start > 1, the burn in should be longer
  heidel <- coda::heidel.diag(x)
  heidel_tidy <- heidel %>%
    purrr::map(unclass) %>%
    purrr::map(tibble::as_tibble, rownames = "var") %>%
    purrr::imap(~ dplyr::mutate(.x, chain = .y, .before = 1)) %>%
    dplyr::bind_rows() %>%
    dplyr::rename_with(~ stringr::str_c(.x, "_heidel"), c(-var, -chain))
  
  # Raftery and Lewis diagnostic
  # Burn-in  Total Lower bound  Dependence
  raftery <- coda::raftery.diag(x)
  raftery_tidy <- raftery %>%
    purrr::map(unclass) %>%
    purrr::modify("resmatrix") %>%
    purrr::map(tibble::as_tibble, rownames = "var") %>%
    purrr::imap(~ dplyr::mutate(.x, chain = .y, .before = 1)) %>%
    dplyr::bind_rows() %>%
    dplyr::rename_with(~ stringr::str_c(.x, "_raftery"), c(-var, -chain))
  
  dplyr::full_join(gelman_tidy, heidel_tidy, raftery_tidy, by = "var")
}

tidy.mcmc <- function(x, vars = "all", ...) {
  
  subset <- subset.mcmc(x, vars)
  
  subset %>%
    purrr::map(tibble::as_tibble) %>%
    purrr::imap(~ dplyr::mutate(.x, chain = .y, .before = 1)) %>%
    dplyr::bind_rows()
}

summary.mcmc <- function(x, vars = "all", bci = FALSE, bci.level = 0.95, ...) {
  
  subset <- subset.mcmc(x, vars)
  
  result <- subset %>%
    summary() %>%
    purrr::pluck("statistics") %>%
    tibble::as_tibble(rownames = "term") %>%
    dplyr::select(term, mean = Mean, sd = SD)
  
  if (bci) {
    bci.min <- (1 - bci.level) / 2
    bci.max <- 1 - bci.min
    tidy <- subset %>% tidy.mcmc() %>% dplyr::select(-chain)
    bci <- tidy %>%
      dplyr::summarize(
        dplyr::across(.fns = ~ quantile(.x, c(bci.min, bci.max))),
        prob = c(bci.min, bci.max)
      ) %>%
      tidyr::pivot_longer(-prob, names_to = "term") %>% 
      tidyr::pivot_wider(names_from = prob, values_from = value)
    result <- dplyr::left_join(result, bci, by = "term")
  }
  
  result
}

summary.jags <- function(x, var, .f, ...) {
  
  result <- x %>% 
    purrr::pluck(var) %>%
    summary(FUN = .f, ...)
  
  result %>%
    purrr::pluck("stat") %>%
    t() %>%
    tibble::as_tibble()
}

generate_seeds <- function(meta_seed, n, range = c(0, 99999)) {
  
  set.seed(meta_seed)
  sample.int(range[2], size = n)
}


# Load reference files
source(file.path(settings$dir, "Flux/tools/reference/seeds.R"))
seeds <- purrr::pluck(
  seeds, settings$site, "patch_model", rlang::as_label(control$flux)
)
source(file.path(settings$dir, "Flux/tools/reference/site_metadata.R"))
md <- purrr::pluck(site_metadata, settings$site)

# Load functions
path_funs <- file.path(settings$dir, "Flux/tools/engine/functions")
source(file.path(path_funs, "clean.R"))
source(file.path(path_funs, "latest_version.R"))
source(file.path(path_funs, "utilities.R"))

# Set working directory and sub-directory structure
wd <- file.path(settings$dir, "Flux", settings$site)

paths <- list(
  flux = file.path(wd, settings$years, "processing", "06_flux_qc", "data"),
  fetch = file.path(wd, settings$years, "processing", "05_footprint", "fetch"),
  cover = file.path(
    wd, settings$years, "processing", "07_footprint_cover", "data"
  ),
  out = file.path(wd, "analysis", "spatial_heterogeneity")
)

tag <- create_tag(settings$site, date = settings$date)

files <- list(
  # Inputs
  flux = purrr::map(paths$flux, ~ latest_version(.x, "flux_qc_full")),
  fetch = purrr::map(paths$fetch, ~ latest_version(.x, "fetch_K15")),
  cover = purrr::map(paths$cover, ~ latest_version(.x, "footprint_cover")),
  chamber = file.path(
    settings$dir, "Gas", c("chamber_ch4_flux.csv", "chamber_co2_flux.csv")
  )
)
# Add output files
files <- append(
  files,
  list() %>%
    append(list(c("model", "samples", "preds"))) %>%
    append(list(c("patch"))) %>%
    append(list(c("inform", "vague"))) %>%
    append(list(tag)) %>%
    expand.grid() %>%
    apply(1, paste, collapse = "_") %>%
    paste0(c(".txt", ".csv", ".csv")) %>%
    paste0(control$flux, "_", .) %>%
    file.path(paths$out, .) %>%
    as.list() %>%
    rlang::set_names(paste0(
      c("model_", "samples_", "preds_"), 
      rep(1, each = 3), rep(c("a", "b"), each = length(rep(1, each = 3)))
    ))
)

library(rjags)
#library(coda)

# Read data
data_flux <- files$flux %>% 
  purrr::map(
    readr::read_csv, guess_max = 7000, 
    col_types = readr::cols(.default = readr::col_guess()), progress = FALSE
  ) %>%
  dplyr::bind_rows()

data_fetch <- files$fetch %>% 
  purrr::map(
    readr::read_csv, guess_max = 7000, 
    col_types = readr::cols(.default = readr::col_guess()), progress = FALSE
  ) %>%
  dplyr::bind_rows()

data_cover <- files$cover %>% 
  purrr::map(
    readr::read_csv, guess_max = 7000, 
    col_types = readr::cols(.default = readr::col_guess()), progress = FALSE
  ) %>%
  dplyr::bind_rows()

data_filter <- data_flux %>%
  dplyr::mutate(
    zm = md$tower_height - md$displacement,
    z0 = calc_z0(ws, ustar, mo_length, zm)
  ) %>%
  #dplyr::full_join(data_fetch, by = "timestamp") %>%
  dplyr::full_join(data_cover, by = "timestamp") %>%
  # Don't need to use phi flag, patch type weights are already filtered
  dplyr::mutate(
    !!control$flux := clean(
      !!control$flux, 
      !!sym(glue::glue("qc_{control$flux}")), 
      !!sym(glue::glue("{control$flux}_unc_flag"))
    )
  ) %>%
  # FILTER: valid/good quality fch4 measurements
  dplyr::filter(!is.na(!!control$flux)) %>%
  # FILTER: daytime
  dplyr::filter(night == 0) %>%
  # FILTER: growing season
  dplyr::filter(dplyr::between(lubridate::month(timestamp), 5, 10)) %>%
  # FILTER: footprint is not concentrated too close to tower
  dplyr::filter(fetch_peak_k15 > 1) %>%
  # FILTER: 85% footprint is within site area 
  dplyr::filter(phi_k15 >= 0.85) %>%
  # FILTER: highest quality sonic anemometer data 
  dplyr::filter(foot_flag == 0) %>%
  # FILTER: sufficient turbulence
  dplyr::filter(ustar > 0.1)
  # FILTER: remove strongly stable conditions
  #dplyr::filter((md$tower_height - md$displacement) / mo_length < 0.4) %>%

# Read chamber data - this is used to specify priors
data_chamber <- files$chamber %>%
  purrr::map(~ readr::read_csv(
    .x, col_types = readr::cols(.default = readr::col_guess())
  )) %>%
  purrr::map(~ dplyr::select(
    .x, site, zone, chamber, ymonth, dplyr::ends_with("flux")
  )) %>% 
  purrr::reduce(
    dplyr::full_join, by = c("site", "zone", "chamber", "ymonth")
  ) %>%
  dplyr::rename(fc = co2_flux, fch4 = ch4_flux) %>%
  dplyr::select(site, zone, flux = !!control$flux) %>%
  tidyr::drop_na() %>%
  dplyr::filter(site == stringr::str_replace(settings$site, "JL", "F")) %>%
  # Convert chamber units (mg C m-2 h-1) to tower units (umol m-2 s-1)
  dplyr::mutate(flux = flux / 12.01 / (60 * 60) * 1000)

# Check distributions
data_chamber %>%
  ggplot(aes(flux)) +
  facet_wrap(~ zone) +
  geom_density() +
  labs(x = rlang::as_label(control$flux))

# Set priors "templates"
priors_a <- data_chamber %>%
  # Recode zone names with parameter names
  dplyr::mutate(
    param = dplyr::recode(zone, A = "beta1", B = "beta2", C = "beta3")
  ) %>%
  dplyr::group_by(param) %>%
  dplyr::summarize(
    mu = mean(flux), sigma = sd(flux), tau = 1 / var(flux), .groups = "drop"
  ) %>%
  as.list() %>%
  # Add vague prior sigma parameters
  append(list(sigma_y = c(control$sigma_prior_min, control$sigma_prior_max)))
priors_b <- priors_a %>%
  purrr::list_modify(mu = rep(0, 3)) %>%
  purrr::list_modify(sigma = rep(100, 3)) %>%
  purrr::list_modify(tau = rep("1.0E-4", 3))


### MODEL 1: PATCH TYPES =======================================================

# Format data for JAGS
# x1, x2, x3 = patch type cover proportions (Dirichlet-distributed)
data_patch <- list(
  y = dplyr::pull(data_filter, !!control$flux), 
  x1 = dplyr::pull(data_filter, wtr), 
  x2 = dplyr::pull(data_filter, veg), 
  x3 = dplyr::pull(data_filter, frs)
)
any(purrr::map_lgl(data_patch, anyNA)) # make sure no missing values

# Model 1A: Patch types (informative priors)

# Prior distributions

# Calculate mu & tau for each patch type
priors_1a <- priors_a
priors_1a

# Multiple chains are used in case results are dependent on initial values
set.seed(seeds$init)
inits_1a <- set_jags_inits(priors_1a, control$n_chains)

# Parameter estimation

# Notes:
# - IMPORTANT: chamber data is NOT normally-distributed - use different 
#   distribution when settings priors??
#   - or apply log transformation? (there are some negative fluxes though)
# - IMPORTANT: what happens when I use random uncertainty as sigma?

# Possible model extensions
# - hierarchical structure with patch types nested within "regions"
#   - i.e. include some form of spatial autocorrelation
# - estimate temperature response parameters for each patch type
# glm <- rstanarm::stan_glm(fch4 ~ 0 + frs + veg + wtr, data = data_filter)
# bayesplot::ppc_dens_overlay(
#   y = glm$y, yrep = rstantools::posterior_predict(glm, draws = 50)
# )
# Set model in JAGS language
model_1a <- set_jags_model(priors_1a, posterior_check = FALSE)
# Write the model to a temporary file for JAGS to access
sink(files$model_1a)
cat(model_1a, fill = TRUE)
sink()

set.seed(seeds$fit)
fit_1a <- rjags::jags.model(
  files$model_1a, data = data_patch, inits = inits_1a, 
  n.chains = control$n_chains, n.adapt = control$n_adapt
)

update(fit_1a, n.iter = control$n_update)

# Sample from the posterior distributions
#rjags::load.module("dic")
samples_1a <- rjags::coda.samples(
  fit_1a, c(paste0("beta", 1:3), "sigma"), n.iter = control$n_iter
)

MCMCvis::MCMCplot(samples_1a, params = c("beta1", "beta2", "beta3"))
summary.mcmc(samples_1a, c("beta1", "beta2", "beta3", "sigma"), bci = TRUE)

#dic_1a <- rjags::dic.samples(fit_1a, n.iter = control$n_iter, type = "pD")
#dic_1a
# Subsample posterior predictions right away to avoid storing massive MCMC list
# set.seed(seeds$pred)
# preds_1a <- samples_1a %>%
#   MCMCvis::MCMCchains(params = "y.sim") %>%
#   magrittr::extract(sample.int(nrow(.), size = 50), ) %>%
#   t() %>% 
#   tibble::as_tibble(.name_repair = "minimal") %>%
#   rlang::set_names(paste0("ysim_", seq(1:ncol(.)))) %>%
#   dplyr::mutate(
#     timestamp = data_filter$timestamp, y = data_patch$y,
#     .before = 1
#   )
# diags_1a <- samples_1a %>%
#   convergence_diags() %>%
#   dplyr::mutate(
#     mean_deviance = sum(as.vector(dic_1a$deviance)),
#     penalty = sum(as.vector(dic_1a$penalty)),
#     penalized_deviance = mean_deviance + penalty
#   )
# Write samples to file
samples_1a %>%
  tidy.mcmc(c("beta1", "beta2", "beta3", "sigma")) %>%
  readr::write_csv(files$samples_1a)
#readr::write_csv(preds_1a, files$preds_1a)

# Credible interval - 95% chance true parameter is within range
#summary(samples_patch)
#MCMCvis::MCMCsummary(samples_patch, round = 4, n.eff = TRUE)

# Plot parameters
samples_1a %>% 
  tidy.mcmc() %>% 
  tidyr::pivot_longer(contains("beta"), values_to = "flux") %>% 
  ggplot(aes(flux, fill = name)) + 
  geom_density(alpha = 0.8) +
  scale_fill_brewer(palette = "YlGnBu", direction = -1)


# Model 1B: Patch types (uninformative priors)

priors_1b <- priors_b
priors_1b

# Initial values are randomly sampled from prior distributions for each chain
set.seed(seeds$init)
inits_1b <- set_jags_inits(priors_1b, control$n_chains)

# Set model in JAGS language
model_1b <- set_jags_model(priors_1b)
# Write the model to a temporary file for JAGS to access
sink(files$model_1b)
cat(model_1b, fill = TRUE)
sink()

set.seed(seeds$fit)
fit_1b <- rjags::jags.model(
  files$model_1b, data = data_patch, inits = inits_1b, 
  n.chains = control$n_chains, n.adapt = control$n_adapt
)

update(fit_1b, n.iter = control$n_update)

# Sample from the posterior distributions
#rjags::load.module("dic")
samples_1b <- rjags::coda.samples(
  fit_1b, c(paste0("beta", 1:3), "sigma"), n.iter = control$n_iter
)
#dic_1b <- rjags::dic.samples(fit_1b, n.iter = control$n_iter, type = "pD")
#dic_1b
# Summarize predictions right away to avoid storing massive MCMC list
# preds_1b <- fit_1b %>%
#   rjags::jags.samples("mu", n.iter = control$n_iter) %>%
#   summary.jags("mu", control$pred_fun, c(0.5, 0.05, 0.95)) %>%
#   dplyr::rename_with(~ paste0("pred_", stringr::str_remove(.x, "%"))) %>%
#   dplyr::mutate(timestamp = data_filter$timestamp, .before = 1) %>%
#   dplyr::bind_cols(tibble::as_tibble(data_patch)) %>%
#   dplyr::mutate(
#     deviance = as.vector(dic_1b$deviance), 
#     penalty = as.vector(dic_1b$penalty)
#   )
# diags_1b <- samples_1b %>%
#   convergence_diags() %>%
#   dplyr::mutate(
#     mean_deviance = sum(as.vector(dic_1b$deviance)),
#     penalty = sum(as.vector(dic_1b$penalty)),
#     penalized_deviance = mean_deviance + penalty
#   )
# Write samples to file
samples_1b %>%
  tidy.mcmc() %>%
  readr::write_csv(files$samples_1b)
#readr::write_csv(preds_1b, files$preds_1b)

# Credible interval - 95% chance true parameter is within range
summary.mcmc(samples_1b, bci = TRUE)

# Plot parameters
samples_1b %>% 
  tidy.mcmc() %>% 
  tidyr::pivot_longer(contains("beta"), values_to = "flux") %>% 
  ggplot(aes(flux, fill = name)) + 
  geom_density(alpha = 0.8) +
  scale_fill_brewer(palette = "YlGnBu", direction = -1)


# Model 1C: Patch types (posterior predictive check)

priors_1c <- priors_a

# Initial values are randomly sampled from prior distributions for each chain
inits_1c <- inits_1a

# Set model in JAGS language
model_1c <- set_jags_model(priors_1c, posterior_check = TRUE)
# Write the model to a temporary file for JAGS to access
sink(file.path(tempdir(), "model_1c.txt"))
cat(model_1c, fill = TRUE)
sink()

set.seed(seeds$fit)
fit_1c <- rjags::jags.model(
  file.path(tempdir(), "model_1c.txt"), data = data_patch, inits = inits_1c, 
  n.chains = control$n_chains, n.adapt = control$n_adapt
)

update(fit_1c, n.iter = control$n_update)

# Sample from the posterior distributions
samples_1c <- rjags::coda.samples(
  fit_1c, c("p.sd", "p.mean", "p.discrep", "p.min", "y.sim"), 
  n.iter = control$n_iter
)
dic_1b <- rjags::dic.samples(fit_1b, n.iter = control$n_iter, type = "pD")
dic_1b
# Summarize predictions right away to avoid storing massive MCMC list
preds_1b <- fit_1b %>%
  rjags::jags.samples("mu", n.iter = control$n_iter) %>%
  summary.jags("mu", control$pred_fun, c(0.5, 0.05, 0.95)) %>%
  dplyr::rename_with(~ paste0("pred_", stringr::str_remove(.x, "%"))) %>%
  dplyr::mutate(timestamp = data_filter$timestamp, .before = 1) %>%
  dplyr::bind_cols(tibble::as_tibble(data_patch)) %>%
  dplyr::mutate(
    deviance = as.vector(dic_1b$deviance), 
    penalty = as.vector(dic_1b$penalty)
  )
diags_1b <- samples_1b %>%
  convergence_diags() %>%
  dplyr::mutate(
    mean_deviance = sum(as.vector(dic_1b$deviance)),
    penalty = sum(as.vector(dic_1b$penalty)),
    penalized_deviance = mean_deviance + penalty
  )
# Write samples to file
samples_1b %>%
  tidy.mcmc() %>%
  readr::write_csv(files$samples_1b)
readr::write_csv(preds_1b, files$preds_1b)

# Credible interval - 95% chance true parameter is within range
summary.mcmc(samples_1b, bci = TRUE)

# Plot parameters
samples_1b %>% 
  tidy.mcmc() %>% 
  tidyr::pivot_longer(contains("beta"), values_to = "flux") %>% 
  ggplot(aes(flux, fill = name)) + 
  geom_density(alpha = 0.8) +
  scale_fill_brewer(palette = "YlGnBu", direction = -1)


### MODEL 2: QUADRANTS =========================================================

# Format data for JAGS
# x1, x2, x3, x4 = quadrant cover proportions (Dirichlet-distributed)
data_quad <- list(
  y = dplyr::pull(data_filter, !!control$flux), 
  x1 = dplyr::pull(data_filter, q1), 
  x2 = dplyr::pull(data_filter, q2), 
  x3 = dplyr::pull(data_filter, q3), 
  x4 = dplyr::pull(data_filter, q4)
)
any(purrr::map_lgl(data_quad, anyNA)) # make sure no missing values

# Prior distributions

# Calculate mu & tau for each quadrant
priors_2b <- list(
  param = paste0("beta", 1:4),
  mu = rep(0, 4),
  sigma = rep(100, 4),
  tau = rep("1.0E-4", 4),
  sigma_y = c(control$sigma_prior_min, control$sigma_prior_max)
)

# Initial values are randomly sampled from prior distributions for each chain
set.seed(seeds$init)
inits_2b <- set_jags_inits(priors_2b, control$n_chains)

# Set model in JAGS language
model_2b <- set_jags_model(priors_2b)
# Write the model to a temporary file for JAGS to access
sink(files$model_2b)
cat(model_2b, fill = TRUE)
sink()

set.seed(seeds$fit)
fit_2b <- rjags::jags.model(
  files$model_2b, data = data_quad, inits = inits_2b, 
  n.chains = control$n_chains, n.adapt = control$n_adapt
)

update(fit_2b, n.iter = control$n_update)

# Sample from the posterior distributions
rjags::load.module("dic")
samples_2b <- rjags::coda.samples(
  fit_2b, c(paste0("beta", 1:4), "sigma", "deviance"), n.iter = control$n_iter
)
#dic_2b <- rjags::dic.samples(fit_2b, n.iter = control$n_iter, type = "pD")
#dic_2b
# Summarize predictions right away to avoid storing massive MCMC list

# diags_2b <- samples_2b %>%
#   convergence_diags() %>%
#   dplyr::mutate(
#     mean_deviance = sum(as.vector(dic_2b$deviance)),
#     penalty = sum(as.vector(dic_2b$penalty)),
#     penalized_deviance = mean_deviance + penalty
#   )
# Write samples to file
samples_2b %>%
  tidy.mcmc(c(paste0("beta", 1:4), "sigma")) %>%
  readr::write_csv(files$samples_2b)
#readr::write_csv(preds_2b, files$preds_2b)

# Credible interval - 95% chance true parameter is within range
summary.mcmc(samples_2b, bci = TRUE)

# Plot parameters
samples_2b %>% 
  tidy.mcmc() %>% 
  tidyr::pivot_longer(contains("beta"), values_to = "flux") %>% 
  ggplot(aes(flux, fill = name)) + 
  geom_density(alpha = 0.8) +
  scale_fill_brewer(palette = "YlGnBu", direction = -1)


### MODEL 3: QUADRANT PATCH TYPES ==============================================

# How are patches distributed across quadrants?
data_filter %>%
  dplyr::select(dplyr::matches("^q\\d+_")) %>%
  tidyr::pivot_longer(
    dplyr::everything(), names_to = c("quad", "class"), names_sep = "_"
  ) %>%
  ggplot2::ggplot(ggplot2::aes(value, fill = class)) +
  ggplot2::facet_wrap(~ quad, scales = "free") +
  ggplot2::geom_density(alpha = 0.8) +
  ggplot2::scale_fill_brewer(palette = "YlGnBu", direction = 1) +
  ggplot2::ylim(0, 50)

