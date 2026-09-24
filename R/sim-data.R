# library(fields)
# library(geoR)
# library(MASS)
# library(spatstat)
# source(here::here("R/bGEV/reparametrized_gev.R"))


#' Simulate preferentially sampled spatial extremes (single period)
#'
#' Simulates a Poisson process on \eqn{[0,1]^2}, a Gaussian field \eqn{S}
#' (called \eqn{Z} in the paper), keeps locations with probability
#' \eqn{\Phi(\beta S(s)/\sigma_S)}, and draws `obs_per_station` GEV maxima
#' per location with median \eqn{\eta + S(s)} and spread \eqn{\exp(\nu)}.
#' Note: responses are simulated from a GEV, not a bGEV.
#'
#' @param beta Preferential-sampling coefficient.
#' @param true_lambda_star Intensity of the underlying Poisson process.
#' @param true_eta,true_nu,true_xi Median, log spread and shape.
#' @param true_sigma2_S,true_rho_S Variance and range of \eqn{S}.
#' @param true_sigma2_T,true_rho_T Variance and range of the random effect on
#'   the spread; used only when `varying_range = TRUE`.
#' @param obs_per_station Number of maxima per location.
#' @param large If `TRUE`, keep at most 150 unobserved locations as
#'   prediction locations; otherwise keep all.
#' @param varying_range Logical; add a random effect to the log spread.
#' @return A list with `y`, `stations`, `obs_coords`, `littlen`, the true
#'   values (`true_q`, `true_log_r`, `true_S_k`, `true_T_k`, `params`) and all
#'   locations (`all_coords`; observed first).
#' @export
make_sim_data_max_stations <-
  function(beta = 2,
           true_lambda_star = 300,
           true_eta = 2,
           true_nu = 1,
           true_sigma2_S = 2,
           true_rho_S = 0.15,
           true_sigma2_T = 0.4, # will be ignored if `varying_range=F`
           true_rho_T = 0.3, # will be ignored if `varying_range=F`
           true_xi = 0.1,
           obs_per_station = 10,
           large = F,
           varying_range = F # decides whether to add random effects to the range parameter
  ) {
    k_true <- c(rpois(1, true_lambda_star))
    true_pp <- spatstat.random::runifpoint(k_true)

    # simulate S_k
    dist.g <- fields::rdist(spatstat.geom::coords(true_pp))
    R_S <-
      exp_cor(dist.g, true_rho_S)
    S_mat <- true_sigma2_S * R_S
    S_k_unordered <- MASS::mvrnorm(1, mu = rep(0, k_true), Sigma = S_mat)

    # Select which coordinates will be kept
    probs <- pnorm(beta * S_k_unordered / sqrt(true_sigma2_S))
    obs_coords <-
      spatstat.geom::coords(true_pp)[(sapply(probs, function(x) {
        rbinom(1, 1, x)
      }) == 1), ]
    idx_keep <- as.integer(rownames(obs_coords))
    S_x <- S_k_unordered[idx_keep]
    n <- dim(obs_coords)[1]

    ## subset_coords
    # idx_pred holds ORIGINAL point indices (previously, with large = TRUE,
    # it held positions within the unobserved subset, so S was taken from the
    # wrong points).
    idx_pred <- which(!(seq_len(k_true) %in% idx_keep))
    if (large == T && length(idx_pred) > 150) {
      idx_pred <- idx_pred[sample(length(idx_pred), 150)]
    }
    pred_coords <- spatstat.geom::coords(true_pp)[idx_pred, ]

    k_pred <- length(idx_pred)

    # some clean up (reorder coordinates and random effects as c(observed, discarded))
    # just to make it easier to keep track of things after.
    S_k <- c(S_x, S_k_unordered[idx_pred])
    # re order coordinates
    all_coords <-
      rbind(obs_coords, pred_coords)

    # Set the median parameter
    q_all <- true_eta + S_k
    q <- q_all[1:n]

    # Simulating the range parameter and the spatial random effects associated with it
    if (varying_range == T) {
      R_T <- exp_cor(fields::rdist(all_coords), range = true_rho_T)
      T_mat <- true_sigma2_T * R_T
      T_k <- MASS::mvrnorm(1, mu = rep(0, nrow(all_coords)), Sigma = T_mat)
      log_r_all <- true_nu + T_k
    } else if (varying_range == F) {
      log_r_all <- true_nu
      T_k <- NULL
    }

    # Shape parameter set in the arguments (true_xi)

    # Simulate data. Keep only selected coordinates
    y_all_mat <- matrix(NA, nrow = (n + k_pred), ncol = obs_per_station)
    for (k in 1:(n + k_pred)) {
      if (varying_range == F) {
        y_all_mat[k, ] <- rgevrep(obs_per_station, q_all[k], exp(log_r_all), true_xi)
      } else if (varying_range == T) {
        y_all_mat[k, ] <- rgevrep(obs_per_station, q_all[k], exp(log_r_all[k]), true_xi)
      }
    }



    y_mat <- y_all_mat[1:n, ]
    y <- as.vector(t(y_mat))
    y_all <- as.vector(t(y_all_mat))
    stations <- rep(1:n, each = obs_per_station)
    stations_all <- rep(1:k, each = obs_per_station)



    return(
      list(
        y = y,
        stations = stations,
        stations_all = stations_all,
        obs_coords = obs_coords,
        littlen = n,
        y_all = y_all,
        true_q = q_all,
        true_log_r = log_r_all,
        true_S_k = S_k,
        true_T_k = T_k,
        all_coords = all_coords,
        params = c(
          "true_lambda_star" = true_lambda_star,
          "true_eta" = true_eta,
          "true_nu" = true_nu,
          "true_sigma2_S" = true_sigma2_S,
          "true_rho_S" = true_rho_S,
          "true_sigma2_T" = true_sigma2_T,
          "true_rho_T" = true_rho_T,
          "true_xi" = true_xi,
          "true_beta" = beta,
          "obs_per_station" = obs_per_station
        )
      )
    )
  }


#' Simulate preferentially sampled space-time extremes
#'
#' As [make_sim_data_max_stations()], with `obs_per_station` periods: each
#' period \eqn{t} has its own \eqn{\eta_t}, \eqn{\nu_t} and spatial effect
#' \eqn{W_t}, added to the shared field \eqn{S}. Calls `set.seed()`
#' internally with a random seed, returned as `gen_seed`.
#'
#' @inheritParams make_sim_data_max_stations
#' @param true_eta_mean,true_nu_mean Means of the per-period \eqn{\eta_t} and
#'   \eqn{\nu_t}.
#' @param true_sigma2_W,true_rho_W Variance and range of \eqn{W_t}.
#' @param obs_per_station Number of periods (one maximum per location per
#'   period).
#' @return As [make_sim_data_max_stations()], plus `years`, `true_W_k` and
#'   `gen_seed`.
#' @export
make_sim_data_max_years <-
  function(
      beta = 2,
      true_lambda_star = 300,
      true_eta_mean = 2,
      true_nu_mean = 1,
      true_sigma2_S = 2,
      true_rho_S = 0.15,
      true_sigma2_W = 0.4,
      true_rho_W = 0.3,
      true_xi = 0.1,
      obs_per_station = 10,
      large = T,
      varying_range = F # decides whether to add random effects to the range parameter
      ) {
    gen_seed <- sample.int(1000, 1)
    set.seed(gen_seed)
    true_eta <- rnorm(obs_per_station, true_eta_mean, 0.5)
    true_nu <- rnorm(obs_per_station, true_nu_mean, 0.2)

    k_true <- c(rpois(1, true_lambda_star))
    true_pp <- spatstat.random::runifpoint(k_true)

    # simulate S_k
    dist.g <- fields::rdist(spatstat.geom::coords(true_pp))
    R_S <-
      exp_cor(dist.g, true_rho_S)
    S_mat <- true_sigma2_S * R_S
    S_k_unordered <- MASS::mvrnorm(1, mu = rep(0, k_true), Sigma = S_mat)

    # Select which coordinates will be kept
    probs <- pnorm(beta * S_k_unordered / sqrt(true_sigma2_S))
    obs_coords <-
      spatstat.geom::coords(true_pp)[
        (sapply(probs, function(x) {
          rbinom(1, 1, x)
        }) ==
          1),
      ]
    idx_keep <- as.integer(rownames(obs_coords))
    S_x <- S_k_unordered[idx_keep]
    n <- dim(obs_coords)[1]

    ## subset_coords
    # idx_pred holds ORIGINAL point indices (previously, with large = TRUE,
    # it held positions within the unobserved subset, so S was taken from the
    # wrong points).
    idx_pred <- which(!(seq_len(k_true) %in% idx_keep))
    if (large == T && length(idx_pred) > 150) {
      idx_pred <- idx_pred[sample(length(idx_pred), 150)]
    }
    pred_coords <- spatstat.geom::coords(true_pp)[idx_pred, ]

    k_pred <- length(idx_pred)

    # some clean up (reorder coordinates and random effects as c(observed, discarded))
    # just to make it easier to keep track of things after.
    S_k <- c(S_x, S_k_unordered[idx_pred])
    # re order coordinates
    all_coords <-
      rbind(obs_coords, pred_coords)

    R_W <-
      exp_cor(dist.g, true_rho_W)
    W_k <- list()
    W_mat <- true_sigma2_W * R_W
    for (j in 1:obs_per_station) {
      W_full <- MASS::mvrnorm(1, mu = rep(0, k_true), Sigma = W_mat)
      # reorder like S_k: observed locations first, then prediction locations
      # (previously left in original point order, misaligned with S_k)
      W_k[[j]] <- W_full[c(idx_keep, idx_pred)]
    }

    q_all <- list()
    q <- list()
    log_r_all <- list()
    for (j in 1:obs_per_station) {
      q_all[[j]] <- true_eta[j] + W_k[[j]] + S_k
      q[[j]] <- q_all[[j]][1:n]
      log_r_all[[j]] <- true_nu[j]
    }

    # Shape parameter set in the arguments (true_xi)

    # Simulate data. Keep only selected coordinates
    y_all_mat <- matrix(NA, nrow = (n + k_pred), ncol = obs_per_station)
    for (j in 1:obs_per_station) {
      for (k in 1:(n + k_pred)) {
        y_all_mat[k, j] <- rgevrep(
          1,
          q_all[[j]][k],
          exp(log_r_all[[j]]),
          true_xi
        )
      }
    }

    y_mat <- y_all_mat[1:n, ]
    y <- as.vector(t(y_mat))
    y_all <- as.vector(t(y_all_mat))
    years <- rep(seq_len(obs_per_station), n) # was rep(1:10, n): only valid for obs_per_station = 10
    stations <- rep(1:n, each = obs_per_station)
    stations_all <- rep(1:k, each = obs_per_station)

    return(
      list(
        y = y,
        stations = stations,
        stations_all = stations_all,
        obs_coords = obs_coords,
        littlen = n,
        years = years,
        y_all = y_all,
        true_q = q_all,
        true_log_r = log_r_all,
        true_S_k = S_k,
        true_W_k = W_k,
        all_coords = all_coords,
        gen_seed = gen_seed,
        params = c(
          "true_lambda_star" = true_lambda_star,
          "true_eta" = true_eta,
          "true_nu" = true_nu,
          "true_sigma2_S" = true_sigma2_S,
          "true_rho_S" = true_rho_S,
          "true_sigma2_W" = true_sigma2_W,
          "true_rho_W" = true_rho_W,
          "true_xi" = true_xi,
          "true_beta" = beta,
          "obs_per_station" = obs_per_station
        )
      )
    )
  }
