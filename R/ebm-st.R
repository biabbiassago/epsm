#' Space-Time Extremes Baseline Model (EBM-ST)
#'
#' Multi-period version of [ebm_mcmc()], and the non-preferential comparison
#' model for [epsm_st_mcmc()]. Each period has its own median, spread and
#' spatial field; the shape and the field's hyperparameters are pooled
#' across periods. There is no point process (no \eqn{\beta}, no
#' \eqn{\lambda^*}).
#'
#' @section Model:
#' For observation \eqn{i} at station \eqn{s_i} in period \eqn{t_i}:
#' \deqn{Y_i \sim \mathrm{bGEV}(\eta_{t_i} + S_{t_i}(s_i),\
#'   \exp(\nu_{t_i}),\ \xi)}
#' with \eqn{S_t(s) = \sigma_S\,\omega_t(s)} and
#' \eqn{\omega_t \mid \rho_S \overset{iid}{\sim} GP(0, R(\cdot;\rho_S))}
#' independently across periods, each defined only at the stations observed
#' in period \eqn{t} (unbalanced panels allowed). \eqn{\xi},
#' \eqn{\sigma^2_S} and \eqn{\rho_S} are shared by all periods.
#'
#' Unlike [epsm_st_mcmc()] there is no time-invariant shared field: the
#' only spatial effect is the period-specific one. For that reason it is
#' **returned under the `W` names** (`W_n`, `sigma2_W`, `rho_W`), so it
#' lines up with the period-specific effect \eqn{W_t} of [epsm_st_mcmc()]
#' when the two fits are compared. Inside the code it is called `S`.
#'
#' @section Priors:
#' \eqn{\eta_t \sim N}, \eqn{\nu_t \sim N}, \eqn{\xi \sim} log-normal,
#' \eqn{\sigma^2_S \sim IG}, \eqn{\rho_S \sim \mathrm{Gamma}}.
#'
#' @section Sampler:
#' Each sweep, for each period \eqn{t}: joint elliptical slice update of
#' \eqn{(\eta_t, \omega_t)}, then adaptive MH for \eqn{\nu_t} (per-period
#' proposal tuning). Then, pooled over periods: adaptive MH for \eqn{\xi}
#' and \eqn{\sigma^2_S}, and MH for \eqn{\rho_S}.
#'
#' @inheritParams ebm_mcmc
#' @param years Vector of length \eqn{N} giving the time period of each
#'   observation. One \eqn{(\eta_t, \nu_t, S_t)} is estimated per unique
#'   value, indexed in order of first appearance. If `NULL`, all data form a
#'   single period (equivalent to [ebm_mcmc()]).
#' @param pred_coords Ignored with a warning: kriging is **not implemented**
#'   for this model. Krige the per-period fields in post-processing from
#'   `W_n`, `sigma2_W` and `rho_W`.
#' @param initial_values Optional named list of starting values; any element
#'   omitted is initialised at random. Recognised names:
#'   \describe{
#'     \item{`eta`, `nu`}{Vectors of length `n_years`.}
#'     \item{`xi`, `sigma2_S`, `rho_S`}{Scalars.}
#'     \item{`S_n`}{List of length `n_years`; element \eqn{t} is a vector of
#'       \eqn{S_t} at the stations observed in period \eqn{t} (natural
#'       scale, ordered as `stations[years == t]`).}
#'   }
#'   If not supplied, each period's \eqn{\nu_t} starts near the log IQR of
#'   that period's data.
#'
#' @return A list of posterior samples (columns / elements index iterations):
#'   \describe{
#'     \item{`eta`, `nu`}{`n_years` x `nsims` matrices.}
#'     \item{`xi`}{Length `nsims`. bGEV shape.}
#'     \item{`W_n`}{Named list (one element per period) of
#'       \eqn{n_t} x `nsims` matrices: \eqn{S_t} at that period's stations.}
#'     \item{`W_pred`}{Always `NULL` (kriging not implemented).}
#'     \item{`sigma2_W`, `rho_W`}{Length `nsims`. Pooled variance and range
#'       of the per-period fields.}
#'     \item{`other_summaries`}{`nu_acc_rate` (`n_years` x `nsims` matrix),
#'       `prev_prop_nu_var` (length `n_years`), `xi_acc_rate`,
#'       `sigma2_W_acc_rate`, `rho_W_acc_rate`, `priors` (a character
#'       summary of the priors used) and `stations`.}
#'   }
#'
#' @seealso [ebm_mcmc()] for the single-period baseline;
#'   [epsm_st_mcmc()] for the preferential-sampling model.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(4649)
#' n <- 30
#' obs_coords <- matrix(runif(2 * n), ncol = 2)
#' stations <- rep(seq_len(n), times = 5)
#' years <- rep(2001:2005, each = n)
#' y <- rnorm(length(stations), 10, 2) # placeholder responses
#' fit <- ebm_st_mcmc(
#'   y = y, obs_coords = obs_coords, stations = stations,
#'   nsims = 5000, years = years, priors = default_priors()
#' )
#' matplot(t(fit$eta), type = "l")
#' }
ebm_st_mcmc <-
  function(y,
           obs_coords,
           stations,
           nsims,
           years = NULL,
           varying_range = F,
           pred_coords = NULL,
           out_file_loc = NULL,
           save = F,
           initial_values = NULL,
           priors = NULL) {
    if (varying_range == TRUE) {
      stop("Not implemented yet")
    }
    if (!length(y) == length(stations)) {
      stop("Stations and y must be two equal length vectors!")
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

    stations_yearly_index <-
      sapply(unique(years), function(x) stations[years == x], simplify = F)
    names(stations_yearly_index) <- names(unique(years))

    sizes_per_year <- sapply(unique(years), function(x) {
      length(stations[years == x])
    })
    names(sizes_per_year) <- unique(years)

    # each observation's position within its own
    # year's omega_j
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
    # initialize chain
    # ------------------------------------------------------------
    eta <- matrix(NA, nrow = n_years, ncol = nsims)
    if (is.null(initial_values$eta)) {
      eta[, 1] <- c(rnorm(n_years, median(y), 0.2))
    } else {
      eta[, 1] <- initial_values$eta
    }

    nu <- matrix(NA, nrow = n_years, ncol = nsims)
    if (is.null(initial_values$nu)) {
      # One starting value PER YEAR, from that year's own data -- see
      # CHANGES-timevar.md section 14 for why a single shared draw
      # recycled across years is worth avoiding.
      overall_iqr <- IQR(y)
      nu[, 1] <- sapply(unique(years), function(yr) {
        iqr_yr <- IQR(y[years == yr])
        if (!is.finite(iqr_yr) || iqr_yr <= 0) {
          iqr_yr <- overall_iqr
        }
        rnorm(1, log(iqr_yr), 0.5)
      })
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
    prev_prop_xi_var <- 0.9

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

    ## omega is from the decomposition: S_j(s) =
    ## sqrt(sigma2_S) * omega_j(s), with omega_j(.) | rho_S ~
    ## GP(0, R_j(.; rho_S)) (unit sill). S_n (below) is the natural-scale
    ## reconstruction, kept for storage/output and for nu/xi's likelihoods.
    dist_obs_coords <- fields::rdist(obs_coords)
    R_S_init <- exp_cor(dist_obs_coords, range = rho_S[1])
    R_S_year_init <-
      sapply(
        unique(idx_years),
        function(x) {
          R_S_init[stations_yearly_index[[x]], stations_yearly_index[[x]]]
        },
        simplify = F
      )

    omega <- vector("list", n_years)
    names(omega) <- unique(years)
    for (yr in seq_len(n_years)) {
      if (!is.null(initial_values$S_n)) {
        omega[[yr]] <- initial_values$S_n[[yr]] / sqrt(sigma2_S[1])
      } else {
        omega[[yr]] <- as.numeric(t(FastGP::rcpp_rmvnorm(
          1,
          S = R_S_year_init[[yr]],
          mu = rep(0, sizes_per_year[yr])
        )))
      }
    }

    S_n <- lapply(sizes_per_year, function(n) {
      matrix(NA, nrow = n, ncol = nsims)
    })
    for (yr in seq_along(S_n)) {
      S_n[[yr]][, 1] <- sqrt(sigma2_S[1]) * omega[[yr]]
    }
    names(S_n) <- unique(years)

    R_S <- list() # rebuilt each sweep, per year, below

    S_pred <- NULL
    if (!is.null(pred_coords)) {
      warning(
        "Kriging for the per-year baseline model is not implemented. Do as post processing."
      )
    }

    #### START RUNNING THE MCMC ITERATIONS #####
    for (i in 2:nsims) {
      if ((i %% 1000) == 0) {
        cat("Iteration:", i, "\n", file = stderr())
      }

      eta_cur_vec <- eta[, i - 1]

      for (j in seq_len(n_years)) {
        y_cur_year <- y[idx_years == j]
        local_stations_cur_year <- stn_idx_in_year[idx_years == j]

        tmp_dist_coords <- fields::rdist(obs_coords[
          stations_yearly_index[[j]],
        ])
        R_S[[j]] <- exp_cor(tmp_dist_coords, range = rho_S[i - 1])
        # Reuse R_S[[j]] just built, rather than recomputing the same
        # distance matrix and correlation function a second time via
        # chol_field_factor().
        chol_R_j <- FastGP::rcppeigen_get_chol(
          R_S[[j]] + diag(1e-6, nrow(R_S[[j]]))
        )

        # Step 1 (per year): joint elliptical slice update of (eta_j,
        # omega_j). Completely independent of every other year's block,
        # given sigma2_S/rho_S -- there is no field shared across years to
        # carry forward here, unlike the preferential timevar model's
        # per-year cycling.
        ess_tmp <-
          sample_eta_S_ess_base(
            eta_cur = eta_cur_vec[j],
            omega_cur = omega[[j]],
            sigma2_S = sigma2_S[i - 1],
            rho_S = rho_S[i - 1],
            nu = nu[j, i - 1],
            xi = xi[i - 1],
            y = y_cur_year,
            stations = local_stations_cur_year,
            coords = obs_coords[stations_yearly_index[[j]], ],
            prior_eta_mean = PRIOR_ETA_MEAN,
            prior_eta_var = PRIOR_ETA_VAR,
            chol_R = chol_R_j
          )
        eta_cur_vec[j] <- ess_tmp$eta
        omega[[j]] <- ess_tmp$omega

        # Natural-scale field for this year, using sigma2_S as it stood
        # BEFORE this sweep's sigma2_S update.
        S_n[[j]][, i] <- sqrt(sigma2_S[i - 1]) * omega[[j]]

        # Scale parameter nu_j (own adaptive-proposal state per year).
        nu_tmp <-
          sample_nu_base(
            nu[j, i - 1],
            eta_cur_vec[j],
            xi[i - 1],
            S_n[[j]][local_stations_cur_year, i],
            y_cur_year,
            prior_nu_mean = PRIOR_NU_MEAN,
            prior_nu_var = PRIOR_NU_VAR,
            iter = i,
            prev_prop_nu_var = prev_prop_nu_var[j],
            full_acc_rate = nu_acc_rate[j, 1:(i - 1)]
          )
        nu[j, i] <- nu_tmp$nu
        nu_acc_rate[j, i] <- nu_tmp$nu_acc
        prev_prop_nu_var[j] <- nu_tmp$prop_nu_var
      }
      eta[, i] <- eta_cur_vec

      # Shape parameter xi: pooled across all years.
      xi_tmp <-
        sample_xi_base_timevar(
          xi[i - 1],
          eta[, i],
          nu[, i],
          lapply(S_n, function(x) x[, i]),
          y,
          flat_idx = flat_idx,
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

      # sigma2_S: pooled across years, holds every omega_j fixed, so only
      # the bGEV likelihood enters.
      sigma2_S_tmp <-
        sample_sigma2_S_base_timevar(
          sigma2_S[i - 1],
          omega,
          eta[, i],
          nu[, i],
          xi[i],
          y,
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

      for (yr in seq_len(n_years)) {
        S_n[[yr]][, i] <- sqrt(sigma2_S[i]) * omega[[yr]]
      }

      # rho_S: pooled across years, holding every omega_j fixed still
      # blocks rho_S from the likelihood entirely. Correctly handles an
      # unbalanced panel (each year's own correlation matrix).
      R_S_inv <- lapply(R_S, function(x) {
        Lx <- FastGP::rcppeigen_get_chol(x + diag(1e-8, nrow(x)))
        chol2inv(t(Lx))
      })
      rho_S_tmp <-
        sample_rho_S_base_timevar(
          rho_S[i - 1],
          R_S,
          R_S_inv,
          omega,
          obs_coords = obs_coords,
          stations_yearly_index = stations_yearly_index,
          rho_prior_a = PRIOR_RHO_S_A,
          rho_prior_b = PRIOR_RHO_S_B
        )
      rho_S[i] <- rho_S_tmp$rho_S
      rho_S_acc_rate[i] <- rho_S_tmp$rho_S_acc
    }

    set_priors <- list(
      eta = paste0("n(", PRIOR_ETA_MEAN, ",", PRIOR_ETA_VAR, ")"),
      nu = paste0("n(", PRIOR_NU_MEAN, ",", PRIOR_NU_VAR, ")"),
      xi = paste0("logn(", PRIOR_XI_MEAN, ",", PRIOR_XI_VAR, ")"),
      sigma2_S = paste0("IG(", PRIOR_SIGMA2_S_A, ",", PRIOR_SIGMA2_S_B, ")"),
      rho_S = paste0("Gamma(", PRIOR_RHO_S_A, ",", PRIOR_RHO_S_B, ")")
    )

    tmp_results <- list(
      eta = eta,
      nu = nu,
      xi = xi,
      W_n = S_n,
      W_pred = S_pred,
      sigma2_W = sigma2_S,
      rho_W = rho_S,
      other_summaries = list(
        nu_acc_rate = nu_acc_rate,
        prev_prop_nu_var = prev_prop_nu_var,
        xi_acc_rate = xi_acc_rate,
        sigma2_W_acc_rate = sigma2_S_acc_rate,
        rho_W_acc_rate = rho_S_acc_rate,
        priors = set_priors,
        stations = stations
      )
    )

    if (save == T) {
      saveRDS(tmp_results, out_file_loc)
    }
    return(tmp_results)
  }
