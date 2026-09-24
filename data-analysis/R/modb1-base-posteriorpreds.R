"This files generates posterior predictives of the median parameters and predictive distribution of y_j(s) for the EB-ST Model"


library(tidyverse)
library(coda)
library(purrr)
source(here::here("data-analysis/R/predictive-utils.R"))

## Note: using a very thinned mcmc output for storage

mcmc_resb1 <- readRDS(here::here("data-analysis/results/pm25_mb1_thinned.rds"))


dat8 <- readRDS(here::here("data-analysis/data/8years_biggrid.rds"))
area_B <- dat8$width_rect * dat8$length_rect

print(dim(dat8$pred_coords))
print(length(unique(dat8$stations)))


# predictions of W

# Create and same S_pred 512 for all these years
W_pred_2024 <- W_pred_by_year_ebm(mcmc_resb1, 2024, dat8$pred_coords, dat8)
# saveRDS(W_pred_2024, here::here("data-analysis/results/W_pred_2024_mt1.rds"))


# W_pred_2023 <- W_pred_by_year(mcmc_rest1, 2023, dat8$pred_coords)
# saveRDS(W_pred_2023, here::here("data-analysis/results/W_pred_2023_mt1.rds"))


# W_pred_2022 <- W_pred_by_year(mcmc_rest1, 2022, dat8$pred_coords)
# saveRDS(W_pred_2022, here::here("data-analysis/results/W_pred_2022_mt1.rds"))


q_pred_2024 <- q_pred_by_year_ebm(mcmc_resb1, 2024, W_pred_2024, dat8$pred_coords)
# saveRDS(q_pred_2024, here::here("data-analysis/results/q_pred_2024_mt1.rds"))

# q_pred_2023 <- q_pred_by_year(mcmc_rest1, 2023, W_pred_2023, dat8$pred_coords)
# saveRDS(q_pred_2023, here::here("data-analysis/results/q_pred_2023_mt1.rds"))

# q_pred_2022 <- q_pred_by_year(mcmc_rest1, 2022, W_pred_2022, dat8$pred_coords)
# saveRDS(q_pred_2022, here::here("data-analysis/results/q_pred_2022_mt1.rds"))


y_pred <- get_y_samples(mcmc_resb1, q_pred_2024, 2024)




#############################################
# ESTIMATED MEDIAN (at observed locations) #
############################################




q_est_2024 <- q_est_by_year_ebm(mcmc_resb1, 2024, dat8$pred_coords, dat8)
saveRDS(q_est_2024, here::here("data-analysis/results/q_est_2024_mt1.rds"))

# q_est_2023 <- q_est_by_year(mcmc_rest1, 2023, dat8$pred_coords,dat8)
# saveRDS(q_est_2023, here::here("data-analysis/results/q_est_2023_mt1.rds"))
#
# q_est_2022 <- q_est_by_year(mcmc_rest1, 2022, dat8$pred_coords,dat8)
# saveRDS(q_est_2022, here::here("data-analysis/results/q_est_2022_mt1.rds"))


y_est_2024 <- get_y_samples(mcmc_resb1, q_est_2024, 2024)
