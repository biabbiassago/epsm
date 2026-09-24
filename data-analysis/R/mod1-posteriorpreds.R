"This files generates posterior predictives of the spatially varying median and reponse y for the EPS-ST Model"

library(tidyverse)
library(coda)
library(purrr)
source(here::here("data-analysis/R/predictive-utils.R"))

## IMPORTANT NOTE: for the preferential model, predictions of the Z(s^*) GP must be made inside the MCMC as the posterior predictive depends on the full augmented point process. Here, the code is to - in post-processing- get predictions of the W_j(s^*) and then of the y_j(s^*) - but the coordinates at which we can predict are only the ones for which posterior predictives are sampled in the mcmc.

## Note 2: using a very thinned mcmc output for storage

# Use
mcmc_rest1 <- readRDS(here::here("data-analysis/results/pm25_m1_thinned.rds"))

dat8 <- readRDS(here::here("data-analysis/data/8years_biggrid.rds"))
area_B <- dat8$width_rect * dat8$length_rect

print(dim(dat8$pred_coords))
print(length(unique(dat8$stations)))


# predictions of W

# Create and same W_pred 512 for all these years
W_pred_2024 <- W_pred_by_year_epsm(mcmc_rest1, 2024, dat8$pred_coords, dat8)
# saveRDS(W_pred_2024, here::here("data-analysis/results/W_pred_2024_mt1.rds"))


## Other years..

# W_pred_2023 <- W_pred_by_year(mcmc_rest1, 2023, dat8$pred_coords,dat8)
# saveRDS(W_pred_2023, here::here("data-analysis/results/W_pred_2023_mt1.rds"))


# W_pred_2022 <- W_pred_by_year(mcmc_rest1, 2022, dat8$pred_coords,dat8)
# saveRDS(W_pred_2022, here::here("data-analysis/results/W_pred_2022_mt1.rds"))



q_pred_2024 <- q_pred_by_year_epsm(mcmc_rest1, 2024, W_pred_2024, dat8$pred_coords)
# saveRDS(q_pred_2024, here::here("data-analysis/results/q_pred_2024_mt1.rds"))

# q_pred_2023 <- q_pred_by_year(mcmc_rest1, 2023, W_pred_2023, dat8$pred_coords)
# saveRDS(q_pred_2023, here::here("data-analysis/results/q_pred_2023_mt1.rds"))

# q_pred_2022 <- q_pred_by_year(mcmc_rest1, 2022, W_pred_2022, dat8$pred_coords)
# saveRDS(q_pred_2022, here::here("data-analysis/results/q_pred_2022_mt1.rds"))


y_pred <- get_y_samples(mcmc_rest1, q_pred_2024, 2024)
library(bayesplot)
mcmc_areas(y_pred[, c(1, 100, 200, 300, 400, 500)]) + ggtitle("Predictive distributions of Y on grid")



#############################################
# ESTIMATED MEDIAN AND Y (at observed locations) #
############################################


q_est_2024 <- q_est_by_year_epsm(mcmc_rest1, 2024, dat8$pred_coords, dat8)
# saveRDS(q_est_2024, here::here("data-analysis/results/q_est_2024_mt1.rds"))

# q_est_2023 <- q_est_by_year(mcmc_rest1, 2023, dat8$pred_coords)
# saveRDS(q_est_2023, here::here("data-analysis/results/q_est_2023_mt1.rds"))
#
# q_est_2022 <- q_est_by_year(mcmc_rest1, 2022, dat8$pred_coords)
# saveRDS(q_est_2022, here::here("data-analysis/results/q_est_2022_mt1.rds"))

y_est_2024 <- get_y_samples(mcmc_rest1, q_est_2024, 2024)
