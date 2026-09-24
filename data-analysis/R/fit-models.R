# This file takes the cleaned dataset for californiam PM2.5 stations
# And runs the EPF-ST and the EB-ST models.

library(coda)
source(here::here("R/epsm-st.R"))
source(here::here("R/samplers-st.R"))
source(here::here("R/ebm-st.R"))
source(here::here("R/utils.R"))


dat8 <- readRDS(here::here("data-analysis/data/8years_biggrid.rds"))
area_B <- dat8$width_rect * dat8$length_rect

# Run Model
NUM_OF_CHAINS <- 1
NSIMS <- 500 ## Change to 100,000 or more iterations. This is just for illustration.

print("Running time varying Preferential model")
# NOTE, for 100,000 iterations this will take roughly 18 hours to run on
# a M3 Macbook
PRIORS_LIST <- make_prior_list()

res_mod1 <- parallel::mclapply(1:NUM_OF_CHAINS, function(x) {
  epsm_st_mcmc(
    y = dat8$y,
    obs_coords = dat8$obs_coords,
    stations = dat8$stations,
    years = dat8$years,
    nsims = NSIMS,
    pred_coords = dat8$pred_coords,
    varying_range = F,
    save = F,
    set_window = owin(xrange = c(0, dat8$width_rect), yrange = c(0, dat8$length_rect)),
    priors = PRIORS_LIST
  )
},
mc.cores = 4
)


params_to_read <- names(res_mod1[[1]][!names(res_mod1[[1]]) %in% c("k", "pred_coords", "other_summaries", "sim_data", "T_n", "sigma2_T", "rho_T")])
res_mod1_mcmc <- make_mcmc_list_obj(res_mod1, params = params_to_read)



print("STARTING BASELINE MODEL NOW.")
res_modb1_timevar <- parallel::mclapply(1:NUM_OF_CHAINS, function(x) {
  ebm_st_mcmc(
    y = dat8$y,
    obs_coords = dat8$obs_coords,
    stations = dat8$stations,
    years = dat8$years,
    nsims = NSIMS,
    pred_coords = dat8$pred_coords,
    varying_range = F,
    save = F
  )
},
mc.cores = 4
)

params_to_read <- names(res_modb1[[1]][!names(res_modb1[[1]]) %in% c("k", "pred_coords", "other_summaries", "sim_data", "T_n", "sigma2_T", "rho_T")])
res_modb1_mcmc <- make_mcmc_list_obj(res_modb1, params = params_to_read)
