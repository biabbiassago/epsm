"Functions for posterior andpredictive/ posterior"


library(tidyverse)
library(coda)
library(purrr)
source(here::here("R/bGEV/bGEVcode.R"))
source(here::here("R/bGEV/allbgev.R"))
source(here::here("R/posterior-preds.R"))
source(here::here("R/utils.R"))

W_pred_by_year_ebm <- function(mod, year, pred_coords, dat) {
  # Posterio predictive for W , EB-ST model
  stations_year <- dat$stations[dat$years == as.character(year)]
  obs_coords_year <- dat$obs_coords[stations_year, ]

  all_params <- colnames(mod[[1]])
  rho_W_mcmc <- unlist(mod[, "rho_W"])
  sigma2_W_mcmc <- unlist(mod[, "sigma2_W"])
  W_n_mcmc <- as.matrix(mod[, startsWith(all_params, paste0("W_n_", year))])
  nsims <- length(sigma2_W_mcmc)
  npreds <- dim(pred_coords)[1]
  W_pred <- matrix(NA, nrow = nsims, ncol = npreds)

  for (j in 1:nsims) {
    if ((j %% 1000) == 0) {
      print(j)
    }
    W_pred[j, ] <- get_S_predlocs(
      rho_W_mcmc[j],
      sigma2_W_mcmc[j],
      W_n_mcmc[j, ],
      pred_coords,
      obs_coords_year
    )
  }
  colnames(W_pred) <- paste0("W_pred_", year, ".", 1:npreds)
  W_pred_mcmc <- mcmc(W_pred)
  return(W_pred_mcmc)
}



q_pred_by_year_ebm <- function(mod, year, pred_W, pred_coords) {
  # Posterio predictive for q , EB-ST model
  year_to_idx <- 2017:2024
  idx <- which(year_to_idx == year)
  med <- paste0("eta.", idx)

  q_pred <-
    mcmc(apply(pred_W, 2, function(x) {
      as.vector(as.matrix(mod[, med])) + x
    }))
  colnames(q_pred) <- paste0("q_pred.", 1:dim(pred_coords)[1])
  return(q_pred)
}


get_y_samples <- function(model_res, q_pred, year) {
  # Predictive distribution for y (both models)
  year_to_idx <- 2017:2024
  year_idx <- which(year_to_idx == year)
  nu_cur <- paste0("nu.", year_idx)
  sb <- exp(as.numeric(unlist(model_res[, nu_cur])))
  xi <- as.numeric(unlist(model_res[, "xi"]))
  nsamples <- dim(q_pred)[1]
  k <- dim(q_pred)[2]
  y_post <- matrix(NA, nrow = nsamples, ncol = k)
  for (i in 1:nsamples) {
    if (i %% 1000 == 0) {
      print(i)
    }
    y_post[i, ] <- rbgev2(1, q_pred[i, ], rep(sb[i], k), rep(xi[i], k))
  }
  colnames(y_post) <- paste0("y_pred.", 1:dim(q_pred)[2])
  return(as.mcmc(y_post))
}



q_pred_by_year_epsm <- function(mod, year, pred_W, pred_coords) {
  # Posterio predictive for q , EPS-ST model
  all_params <- colnames(mod[[1]])
  S_p_mat <- mod[, paste0("S_pred.", 1:dim(pred_coords)[1])]
  raneff_n_mat <- do.call(rbind, S_p_mat) + pred_W
  year_to_idx <- 2017:2024
  idx <- which(year_to_idx == year)
  med <- paste0("eta.", idx)

  q_pred_modt1 <-
    mcmc(apply(raneff_n_mat, 2, function(x) {
      as.vector(as.matrix(mod[, med])) + x
    }))
  colnames(q_pred_modt1) <- paste0("q_pred.", 1:dim(pred_coords)[1])
  return(q_pred_modt1)
}


W_pred_by_year_epsm <- function(mod, year, pred_coords, dat) {
  # Posterio predictive for W , EPS-ST model
  stations_year <- dat$stations[dat$years == as.character(year)]
  obs_coords_year <- dat$obs_coords[stations_year, ]

  all_params <- colnames(mod[[1]])
  rho_W_mcmc <- unlist(mod[, "rho_W"])
  sigma2_W_mcmc <- unlist(mod[, "sigma2_W"])
  W_n_mcmc <- as.matrix(mod[, startsWith(all_params, paste0("W_n_", year))])
  nsims <- length(sigma2_W_mcmc)
  npreds <- dim(pred_coords)[1]
  W_pred <- matrix(NA, nrow = nsims, ncol = npreds)

  for (j in 1:nsims) {
    if ((j %% 1000) == 0) {
      print(j)
    }
    W_pred[j, ] <- get_S_predlocs(
      rho_W_mcmc[j],
      sigma2_W_mcmc[j],
      W_n_mcmc[j, ],
      pred_coords,
      obs_coords_year
    )
  }
  colnames(W_pred) <- paste0("W_pred_", year, ".", 1:npreds)
  W_pred_mcmc <- mcmc(W_pred)
  return(W_pred_mcmc)
}


q_est_by_year_ebm <- function(mod, year, pred_coords, dat) {
  # Posterior distributions of q for median for EB-ST model
  stations_year <- dat$stations[dat$years == as.character(year)]
  obs_coords_year <- dat$obs_coords[stations_year, ]

  all_params <- colnames(mod[[1]])

  raneff_n_mat <- mod[, startsWith(all_params, paste0("W_n_", year))]

  year_to_idx <- 2017:2024
  idx <- which(year_to_idx == year)
  med <- paste0("eta.", idx)

  q_n <-
    mcmc(apply(t(raneff_n_mat), 2, function(x) {
      as.vector(as.matrix(mod[, med])) + x
    }))
  colnames(q_n) <- paste0("q_", stations_year)
  return(q_n)
}
q_est_by_year_epsm <- function(mod, year, pred_coords, dat) {
  # Posterior distributions of q for median for EPS-ST model
  stations_year <- dat$stations[dat$years == as.character(year)]
  obs_coords_year <- dat$obs_coords[stations_year, ]

  all_params <- colnames(mod[[1]])
  S_n_mat <- mod[, paste0("S_n.", stations_year)]
  W_n_mat <- mod[, startsWith(all_params, paste0("W_n_", year))]
  raneff_n_mat <- mcmc.list(Map(`+`, S_n_mat, W_n_mat))
  year_to_idx <- 2017:2024
  idx <- which(year_to_idx == year)
  med <- paste0("eta.", idx)

  q_n <-
    mcmc(apply(t(raneff_n_mat), 2, function(x) {
      as.vector(as.matrix(mod[, med])) + x
    }))
  colnames(q_n) <- paste0("q_", stations_year)
  return(q_n)
}
