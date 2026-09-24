# Constants

area_B <- 1
# NSIMS <- 170000
# BURNIN <- 20000
# THIN <- 1
# NUM_OF_CHAINS <- 1

# Priors

# Gamma
PRIOR_LAMBDA_A <- 400
PRIOR_LAMBDA_B <- 1

# Normal
PRIOR_ETA_MEAN <- 0
PRIOR_ETA_VAR <- 10

# Normal
PRIOR_NU_MEAN <- 0
PRIOR_NU_VAR <- 1

# Log-normal
PRIOR_XI_MEAN <- -2.5
PRIOR_XI_VAR <- 1

# Inverse Gamma
PRIOR_SIGMA2_S_A <- 2
PRIOR_SIGMA2_S_B <- 5

PRIOR_SIGMA2_W_A <- 1
PRIOR_SIGMA2_W_B <- 5

# Inverse Gamma
PRIOR_SIGMA2_T_A <- 1
PRIOR_SIGMA2_T_B <- 5

# Gamma
PRIOR_RHO_S_A <- 200
PRIOR_RHO_S_B <- 1000

# Gamma
PRIOR_RHO_W_A <- 100
PRIOR_RHO_W_B <- 1000

PRIOR_RHO_T_A <- 0
PRIOR_RHO_T_B <- 1


# Normal
PRIOR_BETA_MEAN <- 0
PRIOR_BETA_VAR <- 1

exp_cor <- function(distp, range) {
  return(
    (exp(-(distp / range)))
  )
}


#' Default prior hyperparameters
#'
#' Returns the named list of prior hyperparameters (and `area_B`) defined in
#' `constants.R`, in the form expected by the `priors` argument of
#' [epsm_mcmc()], [epsm_st_mcmc()], [ebm_mcmc()] and [ebm_st_mcmc()]. Modify
#' elements to change a prior.
#'
#' @return Named list.
#' @export
#' @examples
#' pri <- default_priors()
#' pri$PRIOR_BETA_VAR <- 4
default_priors <- function() {
  list(
    area_B = area_B,
    PRIOR_LAMBDA_A = PRIOR_LAMBDA_A, PRIOR_LAMBDA_B = PRIOR_LAMBDA_B,
    PRIOR_ETA_MEAN = PRIOR_ETA_MEAN, PRIOR_ETA_VAR = PRIOR_ETA_VAR,
    PRIOR_NU_MEAN = PRIOR_NU_MEAN, PRIOR_NU_VAR = PRIOR_NU_VAR,
    PRIOR_XI_MEAN = PRIOR_XI_MEAN, PRIOR_XI_VAR = PRIOR_XI_VAR,
    PRIOR_SIGMA2_S_A = PRIOR_SIGMA2_S_A, PRIOR_SIGMA2_S_B = PRIOR_SIGMA2_S_B,
    PRIOR_SIGMA2_W_A = PRIOR_SIGMA2_W_A, PRIOR_SIGMA2_W_B = PRIOR_SIGMA2_W_B,
    PRIOR_SIGMA2_T_A = PRIOR_SIGMA2_T_A, PRIOR_SIGMA2_T_B = PRIOR_SIGMA2_T_B,
    PRIOR_RHO_S_A = PRIOR_RHO_S_A, PRIOR_RHO_S_B = PRIOR_RHO_S_B,
    PRIOR_RHO_W_A = PRIOR_RHO_W_A, PRIOR_RHO_W_B = PRIOR_RHO_W_B,
    PRIOR_RHO_T_A = PRIOR_RHO_T_A, PRIOR_RHO_T_B = PRIOR_RHO_T_B,
    PRIOR_BETA_MEAN = PRIOR_BETA_MEAN, PRIOR_BETA_VAR = PRIOR_BETA_VAR
  )
}
