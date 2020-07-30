
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

# Covariance matrices
# | var1   cov12  cov13 |
# | cov21  var2   cov23 |
# | cov31  cov32  var3  |
# - variance on the diagonal
# - covariance on the off diagonal
# - var1 = sigma1^2
# - cov12 = rho * sigma1 * sigma2
# - rho is the correlation coefficient, [-1, 1]


rm(list = ls())
settings <- list(
  name = "Graham Stewart", # user who ran the script
  email = "grahamstewart12@gmail.com", # user's email contact
  site = "JLN", # three letter site code
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
  # Number of iterations used to choose the sampler and assure optimum mixing of 
  #   the MCMC chain
  n_adapt = 1000,
  # Number of iterations discarded to allow the chain to converge before 
  #   iterations are stored
  n_update = 10000,
  # Number of iterations stored in the chain as samples from the posterior 
  #   distribution
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
    append(list(c("vague", "inform"))) %>%
    append(list(c("nopool", "hier"))) %>%
    append(list(tag)) %>%
    expand.grid() %>%
    apply(1, paste, collapse = "_") %>%
    paste0(c(".txt", ".csv", ".csv")) %>%
    paste0(control$flux, "_", .) %>%
    file.path(paths$out, .) %>%
    as.list() %>%
    rlang::set_names(paste0(
      c("mod", "samps", "preds"),
      rep(1:2, each = 6), rep(c("a", "b"), each = 3)
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
  dplyr::filter(ustar > 0.1)

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

# Set priors "templates"
priors <- data_chamber %>%
  # Recode zone names with parameter names
  dplyr::mutate(
    param = dplyr::recode(zone, A = "beta1", B = "beta2", C = "beta3")
  ) %>%
  dplyr::group_by(param) %>%
  dplyr::summarize(
    mu = mean(flux), sigma = sd(flux), sigma2 = var(flux), tau = 1 / sigma2, 
    .groups = "drop"
  ) %>%
  as.list() %>%
  # Add vague prior sigma parameters
  append(list(sigma_y = c(control$sigma_prior_min, control$sigma_prior_max)))
priors_b <- priors %>%
  purrr::list_modify(mu = rep(0, 3)) %>%
  purrr::list_modify(sigma = rep(100, 3)) %>%
  purrr::list_modify(tau = rep("1.0E-4", 3))

data_grouppatch <- data_filter %>%
  dplyr::select(timestamp, y = !!control$flux, dplyr::matches("^q\\d+_")) %>%
  tidyr::pivot_longer(
    c(-timestamp, -y), names_to = c("q", "c"), values_to = "p", names_sep = "_"
  ) %>%
  tidyr::pivot_wider(names_from = c, values_from = p) %>%
  dplyr::mutate(q = as.integer(stringr::str_remove(q, "q")))

n.j <- 4 # 2 groups
data <- list(
  y = data_grouppatch$y,
  q = data_grouppatch$q,
  x1 = data_grouppatch$wtr,
  x2 = data_grouppatch$veg,
  x3 = data_grouppatch$frs
)


### MODEL 1: NO-POOL ===========================================================

set.seed(seeds$init)
inits1 <- list(
  list(
    beta1 = rep(rnorm(1, priors$mu[1], priors$sigma[1]), n.j),
    beta2 = rep(rnorm(1, priors$mu[2], priors$sigma[2]), n.j), 
    beta3 = rep(rnorm(1, priors$mu[3], priors$sigma[3]), n.j),
    sigma = runif(1, 0, 100)
  ),
  list(
    beta1 = rep(rnorm(1, priors$mu[1], priors$sigma[1]), n.j),
    beta2 = rep(rnorm(1, priors$mu[2], priors$sigma[2]), n.j), 
    beta3 = rep(rnorm(1, priors$mu[3], priors$sigma[3]), n.j),
    sigma = runif(1, 0, 100)
  ),
  list(
    beta1 = rep(rnorm(1, priors$mu[1], priors$sigma[1]), n.j),
    beta2 = rep(rnorm(1, priors$mu[2], priors$sigma[2]), n.j), 
    beta3 = rep(rnorm(1, priors$mu[3], priors$sigma[3]), n.j),
    sigma = runif(1, 0, 100)
  )
)

params1 <- c("beta1", "beta2", "beta3", "sigma")


# Model 1A: vague priors

mod1a <- as.character(glue::glue("
model {{

  # priors
  for (j in 1:{n.j}) {{
    beta1[j] ~ dnorm(0, 0.0001)
    beta2[j] ~ dnorm(0, 0.0001)
    beta3[j] ~ dnorm(0, 0.0001)
  }}
  
  sigma ~ dunif(0, 100)
  tau <- 1 / sigma^2

  # likelihood
  for (i in 1:length(y)) {{
    mu[i] <- beta1[q[i]] * x1[i] + beta2[q[i]] * x2[i] + beta3[q[i]] * x3[i]
    y[i] ~ dnorm(mu[i], tau)
  }}

}}
"))

fit1a <- rjags::jags.model(
  textConnection(mod1a), data, n.chains = 3, n.adapt = 1000, inits = inits1
)

update(fit1a, n.iter = 10000)

samps1a <- rjags::coda.samples(fit1a, params1, n.iter = 10000)

MCMCvis::MCMCplot(samps1a, params = c("beta1", "beta2", "beta3"))
summary(samps1a)

coda::gelman.diag(samps1a, multivariate = FALSE) # if PSRF > 1.05, need more iterations
coda::heidel.diag(samps1a) # if start > 1, burn-in should be longer
coda::raftery.diag(samps1a)

samps1a %>% 
  tidy.mcmc() %>%
  readr::write_csv(files$samps1a)


# Model 1B: Informative priors

mod1b <- as.character(glue::glue("
model {{

  # priors
  for (j in 1:{n.j}) {{
    beta1[j] ~ dnorm({priors$mu[1]}, {priors$tau[1]})
    beta2[j] ~ dnorm({priors$mu[2]}, {priors$tau[2]})
    beta3[j] ~ dnorm({priors$mu[3]}, {priors$tau[3]})
  }}
  
  sigma ~ dunif(0, 100)
  tau <- 1 / sigma^2

  # likelihood
  for (i in 1:length(y)) {{
    mu[i] <- beta1[q[i]] * x1[i] + beta2[q[i]] * x2[i] + beta3[q[i]] * x3[i]
    y[i] ~ dnorm(mu[i], tau)
  }}

}}
"))

params1 <- c("beta1", "beta2", "beta3", "sigma")

set.seed(seeds$fit)
fit1b <- rjags::jags.model(
  textConnection(mod1b), data, n.chains = 3, n.adapt = 1000, inits = inits1
)

update(fit1b, n.iter = 10000)

samps1b <- rjags::coda.samples(fit1b, params1, n.iter = 10000)

MCMCvis::MCMCplot(samps1b, params = c("beta1", "beta2", "beta3"))
summary(samps1b)

coda::gelman.diag(samps1b, multivariate = FALSE) # if PSRF > 1.05, need more iterations
coda::heidel.diag(samps1b) # if start > 1, burn-in should be longer
coda::raftery.diag(samps1b)

samps1b %>% 
  tidy.mcmc() %>%
  readr::write_csv(files$samps1b)


### MODEL 2: RANDOM COEFFICIENTS ===============================================

# No intercepts allowed! All flux must be given to a beta!

set.seed(seeds$init)
B <- purrr::map(1:3, ~ c(
  rnorm(1, priors$mu[1], priors$sigma[1]), 
  rnorm(1, priors$mu[2], priors$sigma[2]), 
  rnorm(1, priors$mu[2], priors$sigma[3])
))

inits2 <- list(
  list(
    B = matrix(B[[1]], nrow = n.j, ncol = 3, byrow = TRUE), 
    mu.beta1 = rnorm(1, priors$mu[1], priors$sigma[1]),
    mu.beta2 = rnorm(1, priors$mu[2], priors$sigma[2]), 
    mu.beta3 = rnorm(1, priors$mu[2], priors$sigma[3]),
    sigma = runif(1, 0, 100)
  ),
  list(
    B = matrix(B[[1]], nrow = n.j, ncol = 3, byrow = TRUE), 
    mu.beta1 = rnorm(1, priors$mu[1], priors$sigma[1]),
    mu.beta2 = rnorm(1, priors$mu[2], priors$sigma[2]), 
    mu.beta3 = rnorm(1, priors$mu[2], priors$sigma[3]),
    sigma = runif(1, 0, 100)
  ),
  list(
    B = matrix(B[[1]], nrow = n.j, ncol = 3, byrow = TRUE), 
    mu.beta1 = rnorm(1, priors$mu[1], priors$sigma[1]),
    mu.beta2 = rnorm(1, priors$mu[2], priors$sigma[2]), 
    mu.beta3 = rnorm(1, priors$mu[2], priors$sigma[3]),
    sigma = runif(1, 0, 100)
  )
)
inits2

params2 <- c(
  "mu.beta1", "mu.beta2", "mu.beta3", "beta1", "beta2", "beta3", 
  "sigma", "Sigma"
)


# Model 2A: vague priors

mod2a <- "
data {

  # scale matrix for Wishart prior on precision matrix
  S[1, 1] <- 1
  S[1, 2] <- 0
  S[1, 3] <- 0
  S[2, 1] <- 0
  S[2, 2] <- 1
  S[2, 3] <- 0
  S[3, 1] <- 0
  S[3, 2] <- 0
  S[3, 3] <- 1
  
}

model {

  # priors for within group model
  sigma ~ dunif(0, 100)
  tau <- 1 / sigma^2

  # likelihood
  
  # model for betas (fixed effects)
  for (i in 1:length(y)) {
    mu[i] <- beta1[q[i]] * x1[i] + beta2[q[i]] * x2[i] + beta3[q[i]] * x3[i]
    y[i] ~ dnorm(mu[i], tau)
  }
  
  # model for group coefficients (random effects)
  for (j in 1:4) {
    
    # group level betas
    beta1[j] <- B[j, 1]
    beta2[j] <- B[j, 2]
    beta3[j] <- B[j, 3]
    
    # mean betas
    B.hat[j, 1] <- mu.beta1
    B.hat[j, 2] <- mu.beta2
    B.hat[j, 3] <- mu.beta3
    
    # group level betas are drawn from beta distributions
    B[j, 1:3] ~ dmnorm(B.hat[j, 1:3], Tau)
  }
  
  # priors
  
  # fixed effects
  mu.beta1 ~ dnorm(0, 0.0001)
  mu.beta2 ~ dnorm(0, 0.0001)
  mu.beta3 ~ dnorm(0, 0.0001)
  
  # inverse of covariance matrix (precision matrix)
  Tau[1:3, 1:3] ~ dwish(S[1:3, 1:3], 4)

  # variance-covariance matrix as derived quantity
  Sigma[1:3, 1:3] <- inverse(Tau[1:3, 1:3])
    
}
"

# 10:17
fit2a <- rjags::jags.model(
  textConnection(mod2a), data, n.chains = 3, n.adapt = 5000, inits = inits2
)

update(fit2a, n.iter = 20000)

samps2a <- rjags::coda.samples(fit2a, params2, n.iter = 50000)

MCMCvis::MCMCplot(samps2a, params = c("beta1", "beta2", "beta3"))
summary(samps2a)

coda::gelman.diag(samps2a, multivariate = FALSE) # if PSRF > 1.05, need more iterations
coda::heidel.diag(samps2a) # if start > 1, burn-in should be longer
coda::raftery.diag(samps2a)

samps2a %>% 
  tidy.mcmc() %>%
  readr::write_csv(files$samps2a)


# Generate informative priors on scale matrix
priors_a
S <- eivtools::get_bugs_wishart_scalemat(c(0.082780029, 0.009466109, 0.021273150))
S$bugs.scalemat

mod2b <- "
data {

  # scale matrix for Wishart prior on precision matrix
  S[1, 1] <- 0.1147473
  S[1, 2] <- 0.000000
  S[1, 3] <- 0.000000
  S[2, 1] <- 0.000000
  S[2, 2] <- 0.01312117
  S[2, 3] <- 0.000000
  S[3, 1] <- 0.000000
  S[3, 2] <- 0.000000
  S[3, 3] <- 0.02967155
  
}

model {

  # priors for within group model
  sigma ~ dunif(0, 100)
  tau <- 1 / sigma^2

  # likelihood
  
  # model for betas (fixed effects)
  for (i in 1:length(y)) {
    mu[i] <- beta1[q[i]] * x1[i] + beta2[q[i]] * x2[i] + beta3[q[i]] * x3[i]
    y[i] ~ dnorm(mu[i], tau)
  }
  
  # model for group coefficients (random effects)
  for (j in 1:4) {
    
    # group level betas
    beta1[j] <- B[j, 1]
    beta2[j] <- B[j, 2]
    beta3[j] <- B[j, 3]
    
    # mean betas
    B.hat[j, 1] <- mu.beta1
    B.hat[j, 2] <- mu.beta2
    B.hat[j, 3] <- mu.beta3
    
    # group level betas are drawn from beta distributions
    B[j, 1:3] ~ dmnorm(B.hat[j, 1:3], Tau)
  }
  
  # priors
  
  # fixed effects
  mu.beta1 ~ dnorm(0.43751894, 12.08021)
  mu.beta2 ~ dnorm(0.23038181, 105.64003)
  mu.beta3 ~ dnorm(0.09959173, 47.00761)
  
  # inverse of covariance matrix (precision matrix)
  Tau[1:3, 1:3] ~ dwish(S[1:3, 1:3], 4)

  # variance-covariance matrix as derived quantity
  Sigma[1:3, 1:3] <- inverse(Tau[1:3, 1:3])
    
}
"

# 9:02
fit2b <- rjags::jags.model(
  textConnection(mod2b), data, n.chains = 3, n.adapt = 5000, inits = inits2
)

update(fit2b, n.iter = 20000)

samps2b <- rjags::coda.samples(fit2b, params2, n.iter = 30000)

MCMCvis::MCMCplot(samps2b, params = c("beta1", "beta2", "beta3"))
summary(samps2b)

coda::gelman.diag(samps2b, multivariate = FALSE) # if PSRF > 1.05, need more iterations
coda::heidel.diag(samps2b) # if start > 1, burn-in should be longer
coda::raftery.diag(samps2b)

samps2b %>% 
  tidy.mcmc() %>%
  readr::write_csv(files$samps2b)






