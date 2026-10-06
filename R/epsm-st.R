#' Space-Time Extremes Preferential Sampling Model (EPSM-ST)
#'
#' Extends [epsm_mcmc()] with a period-specific spatial random effect
#' \eqn{W_j(\cdot)} on the bGEV median. The shared field \eqn{S(\cdot)} still
#' drives both the responses and the preferential point process; \eqn{W_j}
#' affects the responses only.
#'
#' @section Model:
#' For observation \eqn{i} at station \eqn{s_i} in period \eqn{k}:
#' \deqn{Y_{i} \sim \mathrm{bGEV}(\eta_{j} + S(s_i) + W_{j}(s_i),\
#'   \exp(\nu_{t_i}),\ \xi)}
#' \itemize{
#'   \item \eqn{S(s) = \sigma_S\,\omega(s)},
#'     \eqn{\omega \mid \rho_S \sim GP(0, R(\cdot;\rho_S))}: shared,
#'     time-invariant field; enters the point process through
#'     \eqn{\beta\,\omega(s)} as in [epsm_mcmc()].
#'   \item \eqn{W_j(s) = \sigma_W\,\omega_{W,j}(s)},
#'     \eqn{\omega_{W,t} \mid \rho_W \overset{iid}{\sim}
#'     GP(0, R(\cdot;\rho_W))} independently across periods, defined only at
#'     the stations observed in period \eqn{t} (unbalanced years allowed).
#'     \eqn{(\sigma^2_W, \rho_W)} are pooled across periods.
#' }
#'
#' **Notation:** the shared latent field called \eqn{S} here (and `S_n`,
#' `S_pred`, `sigma2_S`, `rho_S` in the code and output) is denoted \eqn{Z}
#' in the paper.
#'
#' @section Sampler:
#' TO DO 
#'
#' @inheritParams epsm_mcmc
#' @param years Vector of length \eqn{N} giving the time period of each
#'   observation. Periods are indexed in order of first appearance.
#' @param initial_values Optional named list of starting values. Accepts
#'   every element listed in [epsm_mcmc()], plus:
#'   \describe{
#'     \item{`sigma2_W`}{Scalar, marginal variance of \eqn{W_j}.}
#'     \item{`rho_W`}{Scalar, range of \eqn{W_t}.}
#'     \item{`W_n`}{List of length `n_years`; element \eqn{j} is a vector of
#'       \eqn{W_t} at the stations observed in period \eqn{j} (natural
#'       scale, ordered as `stations[years == t]`).}
#'   }
#' @param priors Named list of prior hyperparameters. Requires every element
#'   listed in [epsm_mcmc()], plus:
#'   \describe{
#'     \item{`PRIOR_SIGMA2_W_A`, `PRIOR_SIGMA2_W_B`}{\eqn{\sigma^2_W}.}
#'     \item{`PRIOR_RHO_W_A`, `PRIOR_RHO_W_B`}{\eqn{\rho_W}.}
#'   }
#'
#' @return A list with every element returned by [epsm_mcmc()], plus:
#'   \describe{
#'     \item{`W_n`}{Named list (one element per period) of
#'       \eqn{n_t} x `nsims` matrices: \eqn{W_t} at that period's stations.}
#'     \item{`sigma2_W`, `rho_W`}{Length `nsims`. Variance and range of
#'       \eqn{W}.}
#'     \item{`T_n`, `sigma2_T`, `rho_T`}{Placeholders for the varying-range
#'       extension; currently always `NULL`.}
#'   }
#'   Differences in `other_summaries`: `nu_acc_rate` is an `n_years` x
#'   `nsims` matrix and `prev_prop_nu_var` a length-`n_years` vector
#'   (per-period adaptation); adds `sigma2_W_acc_rate`, `rho_W_acc_rate`
#'   and `rho_T_acc_rate` (`NULL`).
#'
#' @seealso [epsm_mcmc()] for the spatial-only model.
#' @importFrom spatstat.geom owin
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(4649)
#' sim_dat <- make_sim_data_max_years(
#'   beta = 2,
#'   obs_per_station = 10, # 10 years per station
#'   large = F,
#'   varying_range = F
#' )
#' y <- sim_dat$y
#' obs_coords <- sim_dat$obs_coords
#' stations <- sim_dat$stations
#' years <- sim_dat$years
#' #' fit <- epsm_st_mcmc(
#'   y = y, obs_coords = obs_coords, stations = stations,
#'   nsims = 5000, years = years, priors = default_priors()
#' )
#' plot(fit$sigma2_W, type = "l")
#' }
epsm_st_mcmc <- function(
    y,
    obs_coords,
    stations,
    nsims,
    years = NULL,
    varying_range = F,
    pred_coords = NULL,
    out_file_loc = NULL,
    save = F,
    sim_data = NULL,
    set_window = owin(),
    initial_values = NULL,
    priors = NULL) {
  # Add some checks.
  if (!length(y) == length(stations)) {
    stop("Stations and y must be two equal length vectors!")
  }
  if (!dim(obs_coords)[1] == length(unique(stations))) {
    stop("There must be as many obs coords as unique stations!")
  }
  littlen <- dim(obs_coords)[1]
  if (is.null(years)) {
    message(
      "No years provided. Assuming all obs are from the same year (single eta and nu parameters)"
    )
    n_years <- 1
    years <- rep(1, length(y))
    idx_years <- years
  } else {
    if (!length(y) == length(years)) {
      stop(
        "length of vector years must be same size as y and stations (or be null)!"
      )
    }
    n_years <- length(unique(years))
    idx_years <- match(years, unique(years))
  }

  ## Set priors from list
  if (is.null(priors)) {
    stop("Set priors (see constants.R file)")
  }
  for (name in names(priors)) {
    assign(name, priors[[name]])
  }

  # Define Stations Index by year
  stations_yearly_index <-
    sapply(
      unique(years),
      function(x) {
        stations[years == x]
      },
      simplify = F
    )
  names(stations_yearly_index) <- names(unique(years))

  # Precompute, once
  year_lengths <- sapply(stations_yearly_index, length)
  year_offsets <- c(0, cumsum(year_lengths))
  stn_idx_in_year <- vapply(
    seq_along(idx_years),
    function(i) {
      match(stations[i], stations_yearly_index[[idx_years[i]]])
    },
    integer(1)
  )
  flat_idx <- year_offsets[idx_years] + stn_idx_in_year
  # ------------------------------------------------------------
  # initialize chains -- initial_values$<name>, when supplied, overrides the
  # default random initialisation (same convention
  # as epsm.R).
  # ------------------------------------------------------------
  if (is.null(initial_values$lambda_star)) {
    lambda_star <- c(rpois(1, littlen * 2))
  } else {
    lambda_star <- c(initial_values$lambda_star)
  }
  k_star <- c(rpois(1, lambda_star[1]))
  pp <- spatstat.random::runifpoint(k_star[1], win = set_window)
  all_coords <- rbind(obs_coords, spatstat.geom::coords(pp))
  k <- c(dim(all_coords)[1])

  if (is.null(initial_values$beta)) {
    beta <- c(runif(1, 1, 3))
  } else {
    beta <- c(initial_values$beta)
  }
  beta_acc_rate <- c(1)
  prev_prop_beta_var <- 0.9

  # time varying eta
  eta <- matrix(NA, nrow = n_years, ncol = nsims)
  if (is.null(initial_values$eta)) {
    eta[, 1] <- c(rnorm(n_years, median(y), 0.5))
  } else {
    eta[, 1] <- initial_values$eta
  }

  # time varying nu
  nu <- matrix(NA, nrow = n_years, ncol = nsims)
  if (is.null(initial_values$nu)) {
    nu[, 1] <- c(rnorm(1, log(IQR(y)), 0.5))
  } else {
    nu[, 1] <- initial_values$nu
  }
  nu_acc_rate <- matrix(NA, nrow = n_years, ncol = nsims)
  nu_acc_rate[, 1] <- 1
  prev_prop_nu_var <- rep(0.5, n_years)

  if (is.null(initial_values$xi)) {
    xi <- c(exp(rnorm(1, -1.2, 0.25)))
  } else {
    xi <- c(initial_values$xi)
  }
  xi_acc_rate <- c(1)
  prev_prop_xi_var <- 0.8

  ## GEV Level 3 pars: shared field S (preferential one)
  if (is.null(initial_values$sigma2_S)) {
    sigma2_S <- c(runif(1, 0.1, 4))
  } else {
    sigma2_S <- c(initial_values$sigma2_S)
  }
  sigma2_S_acc_rate <- c(1)
  prev_prop_sigma2_S_var <- 0.25

  if (is.null(initial_values$rho_S)) {
    rho_S <- c(0.2)
  } else {
    rho_S <- c(initial_values$rho_S)
  }
  rho_S_acc_rate <- c(1)
  prev_prop_rho_S_var <- 0.35

  dist_all_coords <- fields::rdist(all_coords)

  ## omega is such that 
  ## S(s) = sqrt(sigma2_S) * omega(s), with omega(.) | rho_S ~
  ## GP(0, R(.; rho_S))). S_k / S_n are reconstructed from omega
  ## whenever the natural-scale field is actually needed.
  R_S <- exp_cor(dist_all_coords, range = 0.5)
  omega_k <-
    as.numeric(t(FastGP::rcpp_rmvnorm(1, S = R_S, mu = rep(0, nrow(R_S)))))

  # If an initial S_n (natural scale, at the observed stations) is supplied,
  # convert it
  if (!is.null(initial_values$S_n)) {
    omega_k[1:littlen] <- initial_values$S_n / sqrt(sigma2_S[1])
  }

  S_n <- matrix(NA, nrow = littlen, ncol = nsims)
  S_n[, 1] <- sqrt(sigma2_S[1]) * omega_k[1:littlen]

  # Random effects W 
  if (is.null(initial_values$sigma2_W)) {
    sigma2_W <- c(runif(1, 0.1, 3))
  } else {
    sigma2_W <- c(initial_values$sigma2_W)
  }
  sigma2_W_acc_rate <- c(1)

  if (is.null(initial_values$rho_W)) {
    rho_W <- c(0.2)
  } else {
    rho_W <- c(initial_values$rho_W)
  }
  rho_W_acc_rate <- c(1)

  dist_obs_coords <- fields::rdist(obs_coords)
  R_W_init <- exp_cor(dist_obs_coords, range = rho_W[1])
  R_W_n_year_init <-
    sapply(
      unique(idx_years),
      function(x) {
        R_W_init[stations_yearly_index[[x]], stations_yearly_index[[x]]]
      },
      simplify = F
    )

  sizes_per_year <- sapply(unique(years), function(x) {
    length(stations[years == x])
  })
  names(sizes_per_year) <- unique(years)

  # omega_W 
  omega_W <- vector("list", length(sizes_per_year))
  names(omega_W) <- unique(years)
  for (yr in seq_along(sizes_per_year)) {
    if (!is.null(initial_values$W_n)) {
      omega_W[[yr]] <- initial_values$W_n[[yr]] / sqrt(sigma2_W[1])
    } else {
      omega_W[[yr]] <- as.numeric(t(FastGP::rcpp_rmvnorm(
        1,
        S = R_W_n_year_init[[yr]],
        mu = rep(0, sizes_per_year[yr])
      )))
    }
  }

  W_n <- lapply(sizes_per_year, function(n) matrix(NA, nrow = n, ncol = nsims))
  for (yr in seq_along(W_n)) {
    W_n[[yr]][, 1] <- sqrt(sigma2_W[1]) * omega_W[[yr]]
  }
  names(W_n) <- unique(years)
  R_W <- list()

  # Random effects on the range not implemented
  if (varying_range == T) {
    stop("Not implemented yet")
  } else {
    T_n <- NULL
    sigma2_T <- NULL
    rho_T <- NULL
    rho_T_acc_rate <- NULL
  }

  # Kriging stuff
  if (!(is.null(pred_coords))) {
    d11 <- precompute_pred_dists(pred_coords)
    k_pred <- dim(pred_coords)[1]
    S_pred <- matrix(NA, nrow = k_pred, ncol = nsims)
    S_k <- sqrt(sigma2_S[1]) * omega_k
    S_pred[, 1] <-
      get_S_predlocs(
        rho_S[1],
        sigma2_S[1],
        S_k,
        pred_coords,
        all_coords,
        dist_11 = d11
      )
  } else {
    S_pred <- NULL
  }

  for (i in 2:nsims) {
    if ((i %% 1000) == 0) {
      cat("Iteration:", i, "\n", file = stderr())
    }

    # Step 1 : Sample Lambda Star
    lambda_star[i] <-
      sample_lambda_star(
        k[i - 1],
        lambda_prior_a = PRIOR_LAMBDA_A,
        lambda_prior_b = PRIOR_LAMBDA_B,
        area_B = 1
      )

    # Step 2: Simulate Discarded Locations. 
    all_coords_prev <- all_coords
    all_coords_tmp <-
      sample_all_coords(
        lambda_star[i],
        beta[i - 1],
        rho_S[i - 1],
        omega_k,
        obs_coords,
        all_coords_prev,
        owin = set_window,
        area_B = area_B
      )
    all_coords <- all_coords_tmp$all_coords
    omega_x_tilde <- all_coords_tmp$omega_x_tilde

    k[i] <- dim(all_coords)[1]
    if (k[i] <= littlen) {
      all_coords <- obs_coords
      k[i] <- littlen
    }

    # Build the field's correlation matrix and factorise it once per iter.
    chol_R_S <- chol_field_factor(all_coords, rho_S[i - 1])

    # Step 3: joint elliptical slice update of (eta_j, omega)
    eta_cur_vec <- eta[, i - 1]
    omega_k_cur <- c(omega_k[1:littlen], omega_x_tilde)
    W_n_prev <- lapply(W_n, function(x) x[, i - 1])

    for (j in seq_len(n_years)) {
      ess_tmp <-
        sample_eta_j_S_ess_timevar(
          eta_full = eta_cur_vec,
          j = j,
          omega_k_cur = omega_k_cur,
          beta = beta[i - 1],
          sigma2_S = sigma2_S[i - 1],
          rho_S = rho_S[i - 1],
          W_n = W_n_prev,
          nu = nu[, i - 1],
          xi = xi[i - 1],
          y = y,
          littlen = littlen,
          all_coords_new = all_coords,
          stations = stations,
          flat_idx = flat_idx,
          idx_years = idx_years,
          prior_eta_mean = PRIOR_ETA_MEAN,
          prior_eta_var = PRIOR_ETA_VAR,
          chol_R = chol_R_S
        )
      eta_cur_vec[j] <- ess_tmp$eta_j
      omega_k_cur <- ess_tmp$omega_k
    }

    eta[, i] <- eta_cur_vec
    omega_k <- omega_k_cur

    # Natural-scale field for this iter remaining steps
    S_k <- sqrt(sigma2_S[i - 1]) * omega_k
    S_n[, i] <- S_k[1:littlen]

    for (j in unique(years)) {
      idx_pars <- which(unique(years) == j)

      y_cur_year <- y[years == j]
      stations_cur_year <- stations_yearly_index[[idx_pars]]
      S_n_cur <- S_n[stations_cur_year, i]

      # Scale parameter nu (eta already updated above)
      nu_tmp <-
        sample_nu_timevar(
          nu[idx_pars, i - 1],
          eta[idx_pars, i],
          xi[i - 1],
          S_n = S_n_cur,
          W_n = W_n[[idx_pars]][, i - 1],
          y = y_cur_year,
          stations = NULL,
          prior_nu_mean = PRIOR_NU_MEAN,
          prior_nu_var = PRIOR_NU_VAR,
          varying_range = varying_range,
          iter = i,
          prev_prop_nu_var = prev_prop_nu_var[idx_pars],
          full_acc_rate = nu_acc_rate[idx_pars, 1:(i - 1)]
        )
      nu[idx_pars, i] <- nu_tmp$nu
      nu_acc_rate[idx_pars, i] <- nu_tmp$nu_acc
      prev_prop_nu_var[idx_pars] <- nu_tmp$prop_nu_var

      tmp_dist_coords <- fields::rdist(obs_coords[
        stations_yearly_index[[idx_pars]],
      ])
      R_W[[idx_pars]] <- exp_cor(tmp_dist_coords, range = rho_W[i - 1])
      
      # Omega_W for the W_n   
      omega_W[[idx_pars]] <-
        sample_omega_W_ess_timevar(
          omega_W[[idx_pars]],
          R_W[[idx_pars]],
          sigma2_W[i - 1],
          eta[idx_pars, i],
          nu[idx_pars, i],
          xi[i - 1],
          y_cur_year,
          S_n_cur
        )
      W_n[[idx_pars]][, i] <- sqrt(sigma2_W[i - 1]) * omega_W[[idx_pars]]
    }

    # Shape parameter Xi
    xi_tmp <-
      sample_xi_timevar(
        xi[i - 1],
        eta[, i],
        nu[, i],
        S_n[, i],
        lapply(W_n, function(x) x[, i]),
        y,
        stations = stations,
        flat_idx = flat_idx,
        idx_years = idx_years,
        prior_xi_mean = PRIOR_XI_MEAN,
        prior_xi_var = PRIOR_XI_VAR,
        varying_range = varying_range,
        iter = i,
        prev_prop_xi_var = prev_prop_xi_var,
        full_acc_rate = xi_acc_rate
      )
    xi[i] <- xi_tmp$xi
    xi_acc_rate[i] <- xi_tmp$xi_acc
    prev_prop_xi_var <- xi_tmp$prop_xi_var

    # Sample sigma2_W, rho_W (year-specific random effect's hyperparameters
    # -- shared across years)
    R_W_inv <- lapply(R_W, function(x) {
      Lx <- FastGP::rcppeigen_get_chol(x + diag(1e-8, nrow(x)))
      chol2inv(t(Lx))
    })

    sigma2_W_tmp <-
      sample_sigma2_W_timevar(
        sigma2_W[i - 1],
        omega_W,
        eta[, i],
        S_n[, i],
        nu[, i],
        xi[i],
        y,
        stations = stations,
        flat_idx = flat_idx,
        idx_years = idx_years,
        sigma2_W_prior_a = PRIOR_SIGMA2_W_A,
        sigma2_W_prior_b = PRIOR_SIGMA2_W_B
      )
    sigma2_W[i] <- sigma2_W_tmp$sigma2_W
    sigma2_W_acc_rate[i] <- sigma2_W_tmp$sigma2_W_acc

    # Keep the stored, natural-scale W_n consistent with the just-updated
    # sigma2_W[i] -- omega_W itself is unchanged by this step.
    for (yr in seq_along(W_n)) {
      W_n[[yr]][, i] <- sqrt(sigma2_W[i]) * omega_W[[yr]]
    }

    # rho_W
    rho_W_tmp <- sample_rho_base_timevar(
      rho_W[i - 1],
      R_W,
      R_W_inv,
      omega_W,
      obs_coords = obs_coords,
      stations_yearly_index = stations_yearly_index,
      rho_prior_a = PRIOR_RHO_W_A,
      rho_prior_b = PRIOR_RHO_W_B
    )
    rho_W[i] <- rho_W_tmp$rho_W
    rho_W_acc_rate[i] <- rho_W_tmp$rho_W_acc

    # Step 5 : Level 3 spatial pars for the shared field S.
    sigma2_S_tmp <-
      sample_sigma2_S_timevar(
        sigma2_S[i - 1],
        omega_k,
        eta[, i],
        lapply(W_n, function(x) x[, i]),
        nu[, i],
        xi[i],
        y,
        stations = stations,
        flat_idx = flat_idx,
        idx_years = idx_years,
        sigma2_S_prior_a = PRIOR_SIGMA2_S_A,
        sigma2_S_prior_b = PRIOR_SIGMA2_S_B,
        iter = i,
        prev_prop_sigma2_S_var = prev_prop_sigma2_S_var,
        full_acc_rate = sigma2_S_acc_rate
      )
    sigma2_S[i] <- sigma2_S_tmp$sigma2_S
    sigma2_S_acc_rate[i] <- sigma2_S_tmp$sigma2_S_acc
    prev_prop_sigma2_S_var <- sigma2_S_tmp$prop_sigma2_S_var

    # rho_S: holding omega fixed still blocks rho_S from the likelihood,
    rho_S_tmp <-
      sample_rho_S(
        rho_S[i - 1],
        omega_k,
        all_coords,
        rho_prior_a = PRIOR_RHO_S_A,
        rho_prior_b = PRIOR_RHO_S_B,
        iter = i,
        prev_prop_rho_var = prev_prop_rho_S_var,
        full_acc_rate = rho_S_acc_rate,
        chol_R_cur = chol_R_S
      )
    rho_S[i] <- rho_S_tmp$rho_S
    rho_S_acc_rate[i] <- rho_S_tmp$rho_S_acc
    prev_prop_rho_S_var <- rho_S_tmp$prop_rho_var

    if (varying_range == TRUE) {
      stop("not implemented")
    }

    # Step 6: preferential parameter beta.
    beta_tmp <-
      sample_beta(
        beta[i - 1],
        omega_k,
        littlen,
        k[i],
        beta_prior_mean = PRIOR_BETA_MEAN,
        beta_prior_var = PRIOR_BETA_VAR,
        iter = i,
        prev_prop_beta_var = prev_prop_beta_var,
        full_acc_rate = beta_acc_rate
      )
    beta[i] <- beta_tmp$beta
    beta_acc_rate[i] <- beta_tmp$beta_acc
    prev_prop_beta_var <- beta_tmp$prop_beta_var

    # Step 7 : Kriging. Reconstruct S at the fully updated (sigma2_S[i],
    # omega_k) for this sweep.
    if (!(is.null(pred_coords))) {
      S_k <- sqrt(sigma2_S[i]) * omega_k
      S_pred[, i] <-
        get_S_predlocs(
          rho_S[i],
          sigma2_S[i],
          S_k,
          pred_coords,
          all_coords,
          dist_11 = d11
        )
    }
  }

  # House keeping : Set priors and save results
  set_priors <- make_prior_list()
  tmp_results <- list(
    lambda_star = lambda_star,
    k = k,
    eta = eta,
    nu = nu,
    xi = xi,
    beta = beta,
    S_n = S_n,
    W_n = W_n,
    T_n = T_n,
    sigma2_S = sigma2_S,
    sigma2_W = sigma2_W,
    sigma2_T = sigma2_T,
    rho_S = rho_S,
    rho_W = rho_W,
    rho_T = rho_T,
    S_pred = S_pred,
    pred_coords = pred_coords,
    other_summaries = list(
      nu_acc_rate = nu_acc_rate,
      prev_prop_nu_var = prev_prop_nu_var,
      xi_acc_rate = xi_acc_rate,
      beta_acc_rate = beta_acc_rate,
      sigma2_S_acc_rate = sigma2_S_acc_rate,
      sigma2_W_acc_rate = sigma2_W_acc_rate,
      rho_S_acc_rate = rho_S_acc_rate,
      rho_W_acc_rate = rho_W_acc_rate,
      rho_T_acc_rate = rho_T_acc_rate,
      priors = set_priors,
      stations = stations
    )
  )
  if (save == T) {
    saveRDS(tmp_results, out_file_loc)
  }
  return(tmp_results)
}
