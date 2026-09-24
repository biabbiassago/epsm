"Minimal Examples to Run both the EPF Model and the EB Models"
library(epsm)
sim_dat <- make_sim_data_max_stations(beta = 2, obs_per_station = 1)

y <- sim_dat$y
obs_coords <- sim_dat$obs_coords
stations <- sim_dat$stations
# Generate predictions at the "leftover coordinates" from the HPP. This can be any set of coordinates to
# perform prediction at.
pred_coords <- sim_dat$all_coords[(sim_dat$littlen + 1):dim(sim_dat$all_coords)[1], ]


######################################
# ------SPATIAL MODELS -----------#
######################################




# 1. Extremal Preferential Sampling Model

# Create priors. Can use constant.R file.
NSIMS <- 1000
BURNIN <- 10
THIN <- 2
## EPF Model
# Single Chain
mod1 <- epsm_mcmc(
  y = y,
  obs_coords = obs_coords,
  pred_coords = pred_coords,
  stations = stations,
  nsims = NSIMS, # note that longer runs needed for convergence
  priors = PRIOR_LIST
)

## Function to put into coda object

params_to_read <- names(mod1[!names(mod1) %in% c("k", "pred_coords", "other_summaries", "sim_data", "T_n", "sigma2_T", "rho_T")])
# this is then a usual coda mcmc.list object.
mod1_mcmc <- window(make_mcmc_list_obj(list(mod1), burnin = BURNIN, params = params_to_read), thin = THIN)

# Can then treat it as a regular coda mcmc object
summary(mod1_mcmc[, "beta"])
plot(mod1_mcmc[, "eta.1"])

## Note: for the preferential models, predictions are performed in the MCMC as it requires the full augmented field. S_pred are the predictions at given points.
summary(mod1_mcmc[, "S_pred.1"])



###################################################################################################################

# 2. Extremal Baseline Model
# Single Chain
modb1 <- ebm_mcmc(
  y = y,
  obs_coords = obs_coords,
  stations = stations,
  pred_coords = pred_coords,
  nsims = NSIMS,
  priors = default_priors()
)


## If would like to put it into coda object

params_to_read <- names(modb1[!names(modb1) %in% c("k", "pred_coords", "other_summaries", "sim_data", "T_n", "sigma2_T", "rho_T")])
# this is then a usual coda mcmc.list object.
modb1_mcmc <- window(make_mcmc_list_obj(list(modb1), burnin = BURNIN, params = params_to_read), thin = THIN)

# Can then treat it as a regular coda mcmc object
plot(modb1_mcmc[, "eta.1"])

## Note: For spatial only baseline model, predictions can be performed in or out(as post hoc) of the mcmc iterations.
summary(modb1_mcmc[, "S_n.1"])





######################################
# ------SPATIO-TEMPORAL MODELS -----------#
######################################


# Create priors. Can use constant.R file.

# 1. Extremal Preferential Sampling Model - SPATIO TEMPORAL


sim_dat <- make_sim_data_max_years(
  beta = 2,
  obs_per_station = 10, # 10 years per station
  large = F,
  varying_range = F
)

y <- sim_dat$y
obs_coords <- sim_dat$obs_coords
stations <- sim_dat$stations
years <- sim_dat$years
# Generate predictions at the "leftover coordinates" from the HPP. This can be any set of coordinates to
# perform prediction at.
pred_coords <- sim_dat$all_coords[(sim_dat$littlen + 1):dim(sim_dat$all_coords)[1], ]

mod2 <- epsm_st_mcmc(
  y = y,
  obs_coords = obs_coords,
  stations = stations,
  years = years,
  pred_coords = pred_coords,
  nsims = NSIMS,
  priors = default_priors()
)
params_to_read <- names(mod2[!names(mod2) %in% c("k", "pred_coords", "other_summaries", "sim_data", "T_n", "sigma2_T", "rho_T")])
# this is then a usual coda mcmc.list object.
mod2_mcmc <- window(make_mcmc_list_obj(list(mod2), burnin = BURNIN, params = params_to_read), thin = THIN)


# To get value of W : example location 10, year 7
plot(mod2_mcmc[, "W_n_7.10"])

## A note on predictions:
# - Here, for the preferential model, still predict the shared prefereential latent process S_n. However, prediction is not computed for W_n.year . This has to do as post processing, as especially for multiple year it would become very computationally expensive and it's unnecessary to compute it.
# So this works:
plot(mod2_mcmc[, "S_pred.114"])
# But This does NOT!
# plot(mod2_mcmc[,"W_pred_7.114"])

# 2. Baseline Spatio Temporal
modb2 <- ebm_st_mcmc(
  y = y,
  obs_coords = obs_coords,
  stations = stations,
  years = years,
  pred_coords = pred_coords,
  nsims = NSIMS,
  priors = default_priors()
)
# predictions not done!
params_to_read <- names(modb2[!names(modb2) %in% c("k", "pred_coords", "other_summaries", "sim_data", "T_n", "sigma2_T", "rho_T", "S_pred", "W_pred")])
# this is then a usual coda mcmc.list object.

modb2_mcmc <- window(make_mcmc_list_obj(list(modb2), burnin = BURNIN, params = params_to_read), thin = THIN)
