Fitting the EPS and EB Models
================

This vignette shows how to simulate data and fit the four models in
**epsm**:

|                                     | Spatial       | Space-time       |
|-------------------------------------|---------------|------------------|
| Preferential sampling               | `epsm_mcmc()` | `epsm_st_mcmc()` |
| Baseline (no preferential sampling) | `ebm_mcmc()`  | `ebm_st_mcmc()`  |

Throughout, the latent Gaussian field called `S` in the code and output
is denoted $Z$ in the paper.

``` r
library(epsm)
set.seed(4649)
```

## Settings

Chain length, burn-in and thinning used below. These are for
illustration; longer runs are needed for convergence.

``` r
NSIMS <- 1000
BURNIN <- 100
THIN <- 2
```

Priors are passed as a named list. `default_priors()` returns the
defaults; change individual elements as needed:

``` r
priors <- default_priors()
# e.g. a wider prior on the preferential-sampling parameter:
# priors$PRIOR_BETA_VAR <- 4
str(priors[c("PRIOR_BETA_MEAN", "PRIOR_BETA_VAR", "PRIOR_XI_MEAN")])
#> List of 3
#>  $ PRIOR_BETA_MEAN: num 0
#>  $ PRIOR_BETA_VAR : num 1
#>  $ PRIOR_XI_MEAN  : num -2.5
```

The output of every model is a list of chains. `make_mcmc_list_obj()`
turns it into a standard **coda** `mcmc.list`. This helper keeps only
the sampled parameters, dropping bookkeeping variables for these models:

``` r
to_coda <- function(fit, drop = character(0)) {
  skip <- c("k", "pred_coords", "other_summaries", "sim_data",
            "T_n", "sigma2_T", "rho_T", drop)
  params <- setdiff(names(fit), skip)
  params <- params[!vapply(fit[params], is.null, logical(1))]
  window(make_mcmc_list_obj(list(fit), burnin = BURNIN, params = params),
         thin = THIN)
}
```

# Spatial models

## Simulated data

`make_sim_data_max_stations()` simulates a preferentially sampled set of
stations with one block maximum each (`obs_per_station = 1`).

``` r
sim_dat <- make_sim_data_max_stations(beta = 2, obs_per_station = 1)
 
y <- sim_dat$y
obs_coords <- sim_dat$obs_coords
stations <- sim_dat$stations
 
# Predict at the simulated locations that were NOT sampled. Any matrix of
# coordinates can be used here.
pred_coords <- sim_dat$all_coords[(sim_dat$littlen + 1):nrow(sim_dat$all_coords), ]
 
c(stations = sim_dat$littlen, prediction_locations = nrow(pred_coords))
#>             stations prediction_locations 
#>                  186                  144
```

## 1. Extremal preferential sampling model

``` r
mod1 <- epsm_mcmc(
  y = y,
  obs_coords = obs_coords,
  stations = stations,
  pred_coords = pred_coords,
  nsims = NSIMS,
  priors = priors
)
mod1_mcmc <- to_coda(mod1)
 
summary(mod1_mcmc[, "beta"])
plot(mod1_mcmc[, "eta.1"])
```

For the preferential models, predictions are made **inside** the MCMC,
because they need the full augmented field (observed plus thinned
locations). `S_pred.j` is the field at the `j`-th row of `pred_coords`:

``` r
summary(mod1_mcmc[, "S_pred.1"])
```

## 2. Extremal baseline model

Same data, without the point process for the station locations:

``` r
modb1 <- ebm_mcmc(
  y = y,
  obs_coords = obs_coords,
  stations = stations,
  pred_coords = pred_coords,
  nsims = NSIMS,
  priors = priors
)
modb1_mcmc <- to_coda(modb1)
 
plot(modb1_mcmc[, "eta.1"])
summary(modb1_mcmc[, "S_n.1"])
```

For the spatial baseline model, predictions can be made inside the MCMC
(as here, via `pred_coords`) or afterwards from the posterior samples of
`S_n`, `sigma2_S` and `rho_S`.

# Space-time models

## Simulated data

`make_sim_data_max_years()` simulates one maximum per station per year,
with a shared field $S$ plus a year-specific field $W_t$:

``` r
sim_dat <- make_sim_data_max_years(
  beta = 2,
  obs_per_station = 10, # 10 years per station
  large = FALSE,
  varying_range = FALSE
)
 
y <- sim_dat$y
obs_coords <- sim_dat$obs_coords
stations <- sim_dat$stations
years <- sim_dat$years
pred_coords <- sim_dat$all_coords[(sim_dat$littlen + 1):nrow(sim_dat$all_coords), ]
 
table(years)[1:3]
#> years
#>   1   2   3 
#> 190 190 190
```

## 3. Extremal preferential sampling model, space-time

``` r
mod2 <- epsm_st_mcmc(
  y = y,
  obs_coords = obs_coords,
  stations = stations,
  years = years,
  pred_coords = pred_coords,
  nsims = NSIMS,
  priors = priors
)
mod2_mcmc <- to_coda(mod2)
```

Year-specific effects are named `W_n_<year>.<station>`. For example, $W$
in year 7 at the 10th station observed that year:

``` r
plot(mod2_mcmc[, "W_n_7.10"])
```

**Predictions.** The shared field $S$ is predicted inside the MCMC, as
in the spatial model. The year-specific fields $W_t$ are **not**:
kriging every year at every iteration is expensive and unnecessary, so
do it in post-processing from `W_n`, `sigma2_W` and `rho_W`.

``` r
plot(mod2_mcmc[, "S_pred.10"])   # works
# mod2_mcmc[, "W_pred_7.10"]     # does not exist
```

## 4. Extremal baseline model, space-time

The space-time baseline has no shared field: its only spatial effect is
year-specific, and it is returned under the `W` names (`W_n`,
`sigma2_W`, `rho_W`) so it lines up with `epsm_st_mcmc()`. It does no
kriging, so `pred_coords` is ignored (with a warning) and is omitted
here.

``` r
modb2 <- ebm_st_mcmc(
  y = y,
  obs_coords = obs_coords,
  stations = stations,
  years = years,
  nsims = NSIMS,
  priors = priors
)
modb2_mcmc <- to_coda(modb2)
 
summary(modb2_mcmc[, c("sigma2_W", "rho_W")])
```

## Multiple chains

`make_mcmc_list_obj()` takes a list of fits, so several chains can be
combined for convergence diagnostics:

``` r
fits <- lapply(1:3, function(i) {
  ebm_mcmc(y = sim_dat$y[years == 1], obs_coords = obs_coords,
           stations = stations[years == 1], nsims = NSIMS, priors = priors)
})
chains <- make_mcmc_list_obj(fits, burnin = BURNIN, params = c("eta", "xi", "rho_S"))
coda::gelman.diag(chains)
```
