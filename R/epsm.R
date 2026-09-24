#' Extremes Preferential Sampling Model (EPSM) for Spatial Extremes
#'
#' Fits a Bayesian hierarchical model for preferentially sampled spatial
#' extremes by MCMC. Maxima (or minima) follow a blended GEV (bGEV)
#' distribution whose median depends on a latent Gaussian process
#' \eqn{S(\cdot)}; the same field drives the point process that generated the
#' monitoring locations, so station placement is informative about
#' \eqn{S}.
#'
#' @section Model:
#' For observation \eqn{i} at station \eqn{s_i} in period \eqn{t_i}:
#' \deqn{Y_i \sim \mathrm{bGEV}(\eta_{t_i} + S(s_i),\ \exp(\nu_{t_i}),\ \xi)}
#' where the bGEV is parameterised by median (location), spread (scale, on
#' the log scale \eqn{\nu}) and shape \eqn{\xi}. The field is partially
#' whitened, \eqn{S(s) = \sigma_S\,\omega(s)} with
#' \eqn{\omega \mid \rho_S \sim GP(0, R(\cdot;\rho_S))} (exponential
#' correlation, unit sill). Observed locations are the retained points of a
#' thinned Poisson process on `set_window` with dominating intensity
#' \eqn{\lambda^*}; the retention probability depends on
#' \eqn{\beta\,\omega(s)}, so \eqn{\beta} measures the strength of
#' preferential sampling (\eqn{\beta = 0}: no preferential sampling).
#'
#' #'
#' **Notation:** the latent field called \eqn{S} here (and `S_n`, `S_pred`,
#' `sigma2_S`, `rho_S` in the code and output) is denoted \eqn{Z} in the
#' main manuscript.
#'
#'
#' @section Sampler:
#' Each iteeration: (1) Gibbs update of \eqn{\lambda^*}; (2) simulate the
#' thinned-out (unobserved) locations; (3) joint elliptical slice update of
#' \eqn{(\eta, \omega)}; (4) adaptive MH for \eqn{\nu} (per period) and
#' \eqn{\xi}; (5) adaptive MH for \eqn{\sigma^2_S} and \eqn{\rho_S};
#' (6) adaptive MH for \eqn{\beta}; (7) optional kriging of \eqn{S} to
#' `pred_coords`.
#'
#' @param y Numeric vector of length \eqn{N}. Responses (block maxima/minima),
#'   one per station-period observation.
#' @param obs_coords \eqn{n \times 2} numeric matrix of observed station
#'   coordinates. Must lie inside `set_window`.
#' @param stations Integer vector of length \eqn{N}. Row index into
#'   `obs_coords` of the station each element of `y` was recorded at. Must
#'   take values in `1:nrow(obs_coords)` and use every station at least once.
#' @param nsims Integer. Number of MCMC iterations (including the initial
#'   state; no burn-in or thinning is applied).
#' @param years Optional vector of length \eqn{N} giving the time period of
#'   each observation. If `NULL`, all observations are treated as one period
#'   (a single \eqn{\eta} and \eqn{\nu}). Otherwise one \eqn{(\eta_t, \nu_t)}
#'   pair is estimated per unique value, in order of first appearance.
#' @param pred_coords Optional \eqn{m \times 2} matrix of prediction
#'   locations. If supplied, \eqn{S} is kriged to these locations at every
#'   iteration (returned as `S_pred`).
#' @param out_file_loc Character. File path for [saveRDS()] output; used
#'   only when `save = TRUE`.
#' @param save Logical. If `TRUE`, the results list is also written to
#'   `out_file_loc`. Default `FALSE`.
#' @param sim_data Currently unused; reserved for passing the true simulated
#'   values in simulation studies.
#' @param set_window A [spatstat.geom::owin()] observation window. Defaults
#'   to the unit square \eqn{[0,1]^2}. Note: the \eqn{\lambda^*} update
#'   currently assumes a window of area 1.
#' @param initial_values Optional named list of starting values; any element
#'   omitted (or `initial_values = NULL`) is initialised at random.
#'   Recognised names:
#'   \describe{
#'     \item{`lambda_star`}{Scalar, Point process intensity.}
#'     \item{`beta`}{Scalar, preferential-sampling coefficient.}
#'     \item{`eta`}{Vector of length `n_years`, bGEV medians.}
#'     \item{`nu`}{Vector of length `n_years`, log bGEV spreads.}
#'     \item{`xi`}{Scalar, bGEV shape.}
#'     \item{`sigma2_S`}{Scalar, marginal variance of \eqn{S}.}
#'     \item{`rho_S`}{Scalar, range of \eqn{S}.}
#'     \item{`S_n`}{Vector of length \eqn{n}, \eqn{S} at the stations
#'       (natural scale).}
#'   }
#' @param priors Named list of prior hyperparameters (see `constants.R` /
#'   [default_priors()]). Required. Elements are assigned into the function
#'   environment, so names must match exactly:
#'   \describe{
#'     \item{`PRIOR_LAMBDA_A`, `PRIOR_LAMBDA_B`}{\eqn{\lambda^*}.}
#'     \item{`PRIOR_ETA_MEAN`, `PRIOR_ETA_VAR`}{\eqn{\eta_t}.}
#'     \item{`PRIOR_NU_MEAN`, `PRIOR_NU_VAR`}{\eqn{\nu_t}.}
#'     \item{`PRIOR_XI_MEAN`, `PRIOR_XI_VAR`}{\eqn{\xi}.}
#'     \item{`PRIOR_SIGMA2_S_A`, `PRIOR_SIGMA2_S_B`}{\eqn{\sigma^2_S}.}
#'     \item{`PRIOR_RHO_S_A`, `PRIOR_RHO_S_B`}{\eqn{\rho_S}.}
#'     \item{`PRIOR_BETA_MEAN`, `PRIOR_BETA_VAR`}{\eqn{\beta}.}
#'     \item{`area_B`}{Area of the observation window, used when simulating
#'       thinned locations.}
#'   }
#' @param varying_range Logical. If `TRUE`, adds spatial random effects on
#'   the range. **Not yet implemented**; `TRUE` raises an error. Default
#'   `FALSE`.
#'
#' @return A list of posterior samples (columns / elements index iterations):
#'   \describe{
#'     \item{`lambda_star`}{Length `nsims`. Point process intensity.}
#'     \item{`k`}{Length `nsims`. Total points (observed + thinned).}
#'     \item{`eta`}{`n_years` x `nsims` matrix. bGEV medians.}
#'     \item{`nu`}{`n_years` x `nsims` matrix. Log bGEV spreads.}
#'     \item{`xi`}{Length `nsims`. bGEV shape.}
#'     \item{`beta`}{Length `nsims`. Preferential-sampling coefficient.}
#'     \item{`S_n`}{\eqn{n} x `nsims` matrix. \eqn{S} at the stations.}
#'     \item{`sigma2_S`, `rho_S`}{Length `nsims`. GP variance and range.}
#'     \item{`S_pred`}{\eqn{m} x `nsims` matrix of kriged \eqn{S}, or `NULL`.}
#'     \item{`pred_coords`}{As supplied.}
#'     \item{`other_summaries`}{List of MH acceptance-rate traces
#'       (`*_acc_rate`), the final `prev_prop_nu_var`, the prior list
#'       (`priors`) and `stations`.}
#'   }
#'
#' @seealso [epsm_st_mcmc()] for the space-time extension with a
#'   period-specific spatial random effect.
#' @importFrom spatstat.geom owin
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 30
#' sim_dat <- make_sim_data_max_stations(beta = 2, obs_per_station = 1)
#' y <- sim_dat$y
#' obs_coords <- sim_dat$obs_coords
#' stations <- sim_dat$stations
#' fit <- epsm_mcmc(
#'   y = y, obs_coords = obs_coords, stations = stations,
#'   nsims = 5000, years = years, priors = default_priors()
#' )
#' plot(fit$beta, type = "l")
#' }
#'
epsm_mcmc <-
  function(y,
           obs_coords,
           stations,
           nsims,
           years = NULL,
           pred_coords = NULL,
           out_file_loc = NULL,
           save = F,
           sim_data = NULL,
           set_window = owin(),
           initial_values = NULL,
           priors = NULL,
           varying_range = F) {
    if (varying_range == TRUE) {
      stop("Not implemented yet")
    }
    # Add some checks.
    if (!length(y) == length(stations)) {
      stop("Stations and y must be two equal length vectors!")
    }
    if (!dim(obs_coords)[1] == length(unique(stations))) {
      stop("There must be as many obs coords as unique stations!")
    }

    ## Set priors from list
    if (is.null(priors)) {
      stop("Set priors (see constants.R file)")
    }
    for (name in names(priors)) {
      assign(name, priors[[name]])
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
      idx_years <- as.numeric(as.factor(years))
    }

    # ------------------------------------------------------------
    # initialize chains
    # ------------------------------------------------------------

    # augmented locations.
    if (is.null(initial_values$lambda_star)) {
      lambda_star <- c(rpois(1, littlen * 2))
    } else {
      lambda_star <- c(initial_values$lambda_star)
    }

    k_star <- c(rpois(1, lambda_star[1]))

    # Generate initial augmented locations
    pp <- spatstat.random::runifpoint(
      k_star[1],
      win = set_window
    )

    all_coords <- rbind(
      obs_coords,
      spatstat.geom::coords(pp)
    )

    k <- c(dim(all_coords)[1])

    if (is.null(initial_values$beta)) {
      beta <- c(runif(1, 1, 3))
    } else {
      beta <- c(initial_values$beta)
    }
    beta_acc_rate <- c(1)
    prev_prop_beta_var <- 0.9

    eta <- matrix(
      NA,
      nrow = n_years,
      ncol = nsims
    )

    if (is.null(initial_values$eta)) {
      eta[, 1] <- c(rnorm(n_years, median(y), 0.2))
    } else {
      eta[, 1] <- initial_values$eta
    }

    eta_acc_rate <- c(1)
    prev_prop_eta_var <- 0.3

    nu <- matrix(
      NA,
      nrow = n_years,
      ncol = nsims
    )
    if (is.null(initial_values$nu)) {
      nu[, 1] <- c(rnorm(1, log(IQR(y)), 0.5))
    } else {
      nu[, 1] <- initial_values$nu
    }

    nu_acc_rate <- c(1)
    prev_prop_nu_var <- 0.5

    if (is.null(initial_values$xi)) {
      xi <- c(exp(rnorm(1, -1.2, 0.25)))
    } else {
      xi <- c(initial_values$xi)
    }

    xi_acc_rate <- c(1)
    prev_prop_xi_var <- 0.8

    ## GEV Level 3 pars

    if (is.null(initial_values$sigma2_S)) {
      sigma2_S <- c(runif(1, 0.1, 3))
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

    ## GEV Random Effects on Median. omega is the decomposed S(s) = sqrt(sigma2_S) * omega(s), with
    ## omega(.) | rho_S ~ GP(0, R(.; rho_S)), unit sill). S_k / S_n are
    ## reconstructed from omega whenever the natural-scale field is
    ## needed (storage, nu/xi, kriging)
    R_S <-
      exp_cor(dist_all_coords, range = 0.5)
    omega_k <-
      as.numeric(t(FastGP::rcpp_rmvnorm(
        1,
        S = R_S,
        mu = rep(0, nrow(R_S))
      )))

    S_n <- matrix(NA, nrow = littlen, ncol = nsims)
    S_k <- matrix(NA, nrow = k, ncol = nsims)
    S_n[, 1] <- sqrt(sigma2_S[1]) * omega_k[1:littlen]
    if (!is.null(initial_values$S_n)) {
      S_k[1:littlen] <- initial_values$S_n
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
    } else if (is.null(pred_coords)) {
      S_pred <- NULL
    }

    #### START RUNNING THE MCMC ITERATIONS #####
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

      # Step 2: Simulate Discarded Locations. Works entirely on the
      # correlation (omega) scale -- sigma2_S never enters the augmentation
      # or thinning step under partial whitening.
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

      # Build the field's correlation matrix and factorise it ONCE per iteration.
      # The same factor is shared by the joint (eta, omega) elliptical slice
      # update and the rho_S update.
      chol_R_S <- chol_field_factor(all_coords, rho_S[i - 1])

      # Step 3: joint elliptical slice update of (eta, omega) at all
      # locations. (Replaces the old field-only ESS step plus the
      # separate eta Metropolis step).
      ess_tmp <-
        sample_eta_S_ess(
          eta_cur = eta[, i - 1],
          omega_n_cur = omega_k[1:littlen],
          omega_x_tilde = omega_x_tilde,
          beta = beta[i - 1],
          sigma2_S = sigma2_S[i - 1],
          rho_S = rho_S[i - 1],
          nu = nu[, i - 1],
          xi = xi[i - 1],
          y = y,
          littlen = littlen,
          all_coords_new = all_coords,
          stations = stations,
          idx_years = idx_years,
          prior_eta_mean = PRIOR_ETA_MEAN,
          prior_eta_var = PRIOR_ETA_VAR,
          chol_R = chol_R_S
        )
      eta[, i] <- ess_tmp$eta
      omega_k <- ess_tmp$omega_k

      # Natural-scale field recompute
      S_k <- sqrt(sigma2_S[i - 1]) * omega_k
      S_n[, i] <- S_k[1:littlen]

      ## Step 4: GEV level 2 parameters (scale nu; eta already updated above)
      for (j in unique(years)) {
        idx_pars <- which(unique(years) == j)

        y_cur_year <- y[years == j]
        stations_cur_year <- stations[years == j]
        S_n_cur <- S_n[stations_cur_year, i]

        nu_tmp <-
          sample_nu(
            nu[idx_pars, i - 1],
            eta[idx_pars, i],
            xi[i - 1],
            S_n_cur,
            y_cur_year,
            prior_nu_mean = PRIOR_NU_MEAN,
            prior_nu_var = PRIOR_NU_VAR,
            iter = i,
            prev_prop_nu_var = prev_prop_nu_var,
            full_acc_rate = nu_acc_rate
          )
        nu[idx_pars, i] <- nu_tmp$nu
        nu_acc_rate[i] <- nu_tmp$nu_acc
        prev_prop_nu_var <- nu_tmp$prop_nu_var
      }

      # Shape parameter Xi
      xi_tmp <-
        sample_xi(
          xi[i - 1],
          eta[, i],
          nu[, i],
          S_n[, i],
          y,
          stations = stations,
          idx_years = idx_years,
          prior_xi_mean = PRIOR_XI_MEAN,
          prior_xi_var = PRIOR_XI_VAR,
          iter = i,
          prev_prop_xi_var = prev_prop_xi_var,
          full_acc_rate = xi_acc_rate
        )
      xi[i] <- xi_tmp$xi
      xi_acc_rate[i] <- xi_tmp$xi_acc
      prev_prop_xi_var <- xi_tmp$prop_xi_var

      # Step 5 : Level 3 spatial pars for the field.
      # sigma2_S: holds omega fixed, so only the bGEV likelihood enters --
      sigma2_S_tmp <-
        sample_sigma2_S(
          sigma2_S[i - 1],
          omega_k,
          eta[, i],
          nu[, i],
          xi[i],
          y,
          stations = stations,
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
      # rho_S
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
      # omega_k).
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

    set_priors <- make_prior_list()
    # save everything to a list.
    tmp_results <- list(
      lambda_star = lambda_star,
      k = k,
      eta = eta,
      nu = nu,
      xi = xi,
      beta = beta,
      S_n = S_n,
      sigma2_S = sigma2_S,
      rho_S = rho_S,
      S_pred = S_pred,
      pred_coords = pred_coords,
      other_summaries = list(
        nu_acc_rate = nu_acc_rate,
        prev_prop_nu_var = prev_prop_nu_var,
        xi_acc_rate = xi_acc_rate,
        beta_acc_rate = beta_acc_rate,
        sigma2_S_acc_rate = sigma2_S_acc_rate,
        rho_S_acc_rate = rho_S_acc_rate,
        priors = set_priors,
        stations = stations
      )
    )
    if (save == T) {
      saveRDS(tmp_results, out_file_loc)
    }
    return(tmp_results)
  }
