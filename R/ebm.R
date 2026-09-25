#' Extremes Baseline Model (EBM) for Spatial Extremes
#'
#' Fits the non-preferential comparison model for [epsm_mcmc()]: a bGEV
#' model for maxima (or minima) with a latent Gaussian process on the
#' median, but **no** point process for the station locations (no
#' \eqn{\beta}, no \eqn{\lambda^*}). Station locations are treated as fixed
#' and uninformative. Single period only: one \eqn{\eta}, one \eqn{\nu}, one
#' field.
#'
#' @section Model:
#' For observation \eqn{i} at station \eqn{s_i}:
#' \deqn{Y_i \sim \mathrm{bGEV}(\eta + S(s_i),\ \exp(\nu),\ \xi)}
#' where the bGEV is parameterised by median (location), spread (scale, on
#' the log scale \eqn{\nu}) and shape \eqn{\xi}. The field is partially
#' whitened, \eqn{S(s) = \sigma_S\,\omega(s)} with
#' \eqn{\omega \mid \rho_S \sim GP(0, R(\cdot;\rho_S))} (exponential
#' correlation, unit sill), defined at the observed stations only.
#'
#' **Notation:** the latent field called \eqn{S} here (and `S_n`, `S_pred`,
#' `sigma2_S`, `rho_S` in the code and output) is denoted \eqn{W} in the
#' paper.
#'
#'
#' @section Sampler:
#' Each iteration: (1) joint elliptical slice update of \eqn{(\eta, \omega)};
#' (2) adaptive MH for \eqn{\nu}, \eqn{\xi}, \eqn{\sigma^2_S} and
#' \eqn{\rho_S}; (3) optional kriging of \eqn{S} to `pred_coords`.
#'
#' @param y Numeric vector of length \eqn{N}. Responses (block maxima/minima).
#' @param obs_coords \eqn{n \times 2} numeric matrix of observed station
#'   coordinates.
#' @param stations Integer vector of length \eqn{N}. Row index into
#'   `obs_coords` of the station each element of `y` was recorded at. Must
#'   take values in `1:nrow(obs_coords)` and use every station at least once.
#' @param nsims Integer. Number of MCMC iterations (including the initial
#'   state; no burn-in or thinning is applied).
#' @param varying_range Logical. If `TRUE`, adds spatial random effects on
#'   the range. **Not yet implemented**; `TRUE` raises an error. Default
#'   `FALSE`.
#' @param pred_coords Optional \eqn{m \times 2} matrix of prediction
#'   locations. If supplied, \eqn{S} is kriged to these locations at every
#'   iteration (returned as `S_pred`).
#' @param out_file_loc Character. File path for [saveRDS()] output; used
#'   only when `save = TRUE`.
#' @param save Logical. If `TRUE`, the results list is also written to
#'   `out_file_loc`. Default `FALSE`.
#' @param initial_values Optional named list of starting values; any element
#'   omitted (or `initial_values = NULL`) is initialised at random.
#'   Recognised names: `eta`, `nu`, `xi`, `sigma2_S`, `rho_S` (scalars) and
#'   `S_n` (vector of length \eqn{n}, \eqn{S} at the stations, natural
#'   scale).
#' @param priors Named list of prior hyperparameters (see `constants.R`).
#'   Required. Elements are assigned into the function environment, so names
#'   must match exactly:
#'   \describe{
#'     \item{`PRIOR_ETA_MEAN`, `PRIOR_ETA_VAR`}{Normal prior on \eqn{\eta}.}
#'     \item{`PRIOR_NU_MEAN`, `PRIOR_NU_VAR`}{Normal prior on \eqn{\nu}.}
#'     \item{`PRIOR_XI_MEAN`, `PRIOR_XI_VAR`}{Log-normal prior on \eqn{\xi}.}
#'     \item{`PRIOR_SIGMA2_S_A`, `PRIOR_SIGMA2_S_B`}{Inverse-gamma prior on
#'       \eqn{\sigma^2_S}.}
#'     \item{`PRIOR_RHO_S_A`, `PRIOR_RHO_S_B`}{Gamma prior on \eqn{\rho_S}.}
#'   }
#'   Other elements (e.g. the preferential-sampling priors) are ignored, so
#'   the same list can be shared with [epsm_mcmc()].
#'
#' @return A list of posterior samples (columns / elements index iterations):
#'   \describe{
#'     \item{`eta`, `nu`}{`1` x `nsims` matrices (same shape as the
#'       single-period output of [epsm_mcmc()]).}
#'     \item{`xi`}{Length `nsims`. bGEV shape.}
#'     \item{`S_n`}{\eqn{n} x `nsims` matrix. \eqn{S} at the stations.}
#'     \item{`S_pred`}{\eqn{m} x `nsims` matrix of kriged \eqn{S}, or `NULL`.}
#'     \item{`sigma2_S`, `rho_S`}{Length `nsims`. GP variance and range.}
#'     \item{`other_summaries`}{List of MH acceptance-rate traces
#'       (`nu_acc_rate`, `xi_acc_rate`, `sigma2_W_acc_rate`,
#'       `rho_W_acc_rate` -- note: these last two are the \eqn{\sigma^2_S}
#'       and \eqn{\rho_S} rates), the final `prev_prop_nu_var`, `priors` (a
#'       character summary of the priors used) and `stations`.}
#'   }
#'
#' @seealso [ebm_st_mcmc()] for the multi-period baseline;
#'   [epsm_mcmc()] for the preferential-sampling model.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(4649)
#' n <- 30
#' obs_coords <- matrix(runif(2 * n), ncol = 2)
#' stations <- seq_len(n)
#' y <- rnorm(n, 10, 2) # placeholder responses
#' fit <- ebm_mcmc(
#'   y = y, obs_coords = obs_coords, stations = stations,
#'   nsims = 5000, priors = default_priors()
#' )
#' plot(fit$rho_S, type = "l")
#' }
ebm_mcmc <-
  function(y,
           obs_coords,
           stations,
           nsims,
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

    # ------------------------------------------------------------
    # initialize chain
    # ------------------------------------------------------------
    if (is.null(initial_values$eta)) {
      eta <- c(rnorm(1, median(y), 0.2))
    } else {
      eta <- c(initial_values$eta)
    }

    if (is.null(initial_values$nu)) {
      nu <- c(rnorm(1, log(IQR(y)), 0.5))
    } else {
      nu <- c(initial_values$nu)
    }
    nu_acc_rate <- c(1)
    prev_prop_nu_var <- 0.5

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
    prev_prop_rho_S_var <- 0.35

    ## omega S(s) = sqrt(sigma2_S)*omega(s),
    ## with omega(.) | rho_S ~ GP(0, R(.; rho_S)), unit sill. S_n (below) is
    ## the natural-scale reconstruction, kept for storage/output and for
    ## nu/xi's likelihoods.
    dist_obs_coords <- fields::rdist(obs_coords)
    R_S_init <- exp_cor(dist_obs_coords, range = rho_S[1])
    omega <-
      as.numeric(t(FastGP::rcpp_rmvnorm(1, S = R_S_init, mu = rep(0, littlen))))

    if (!is.null(initial_values$S_n)) {
      omega <- initial_values$S_n / sqrt(sigma2_S[1])
    }

    S_n <- matrix(NA, nrow = littlen, ncol = nsims)
    S_n[, 1] <- sqrt(sigma2_S[1]) * omega

    # Kriging stuff
    if (!(is.null(pred_coords))) {
      d11 <- precompute_pred_dists(pred_coords)
      k_pred <- dim(pred_coords)[1]
      S_pred <- matrix(NA, nrow = k_pred, ncol = nsims)
      S_pred[, 1] <-
        get_S_predlocs(rho_S[1], sigma2_S[1], S_n[, 1], pred_coords, obs_coords, dist_11 = d11)
    } else {
      S_pred <- NULL
    }

    #### START RUNNING THE MCMC ITERATIONS #####
    for (i in 2:nsims) {
      if ((i %% 10000) == 0) {
        cat("Iteration:", i, "\n", file = stderr())
      }

      # Build the field's correlation matrix and factorise it ONCE per iter.
      # Shared by the joint (eta, omega) elliptical slice update and
      # the rho_S update -- neither needs R^{-1} explicitly.
      chol_R_S <- chol_field_factor(obs_coords, rho_S[i - 1])

      # Step 1: joint elliptical slice update of (eta, omega)
      ess_tmp <-
        sample_eta_S_ess_base(
          eta_cur = eta[i - 1],
          omega_cur = omega,
          sigma2_S = sigma2_S[i - 1],
          rho_S = rho_S[i - 1],
          nu = nu[i - 1],
          xi = xi[i - 1],
          y = y,
          stations = stations,
          coords = obs_coords,
          prior_eta_mean = PRIOR_ETA_MEAN,
          prior_eta_var = PRIOR_ETA_VAR,
          chol_R = chol_R_S
        )
      eta[i] <- ess_tmp$eta
      omega <- ess_tmp$omega

      # Reassign natural scale.
      S_k <- sqrt(sigma2_S[i - 1]) * omega
      S_n[, i] <- S_k

      # Scale parameter nu
      nu_tmp <-
        sample_nu_base(
          nu[i - 1],
          eta[i],
          xi[i - 1],
          S_n[stations, i],
          y,
          prior_nu_mean = PRIOR_NU_MEAN,
          prior_nu_var = PRIOR_NU_VAR,
          iter = i,
          prev_prop_nu_var = prev_prop_nu_var,
          full_acc_rate = nu_acc_rate
        )
      nu[i] <- nu_tmp$nu
      nu_acc_rate[i] <- nu_tmp$nu_acc
      prev_prop_nu_var <- nu_tmp$prop_nu_var

      # Shape parameter Xi
      xi_tmp <-
        sample_xi_base(
          xi[i - 1],
          eta[i],
          nu[i],
          S_n[, i],
          y,
          stations = stations,
          prior_xi_mean = PRIOR_XI_MEAN,
          prior_xi_var = PRIOR_XI_VAR,
          iter = i,
          prev_prop_xi_var = prev_prop_xi_var,
          full_acc_rate = xi_acc_rate
        )
      xi[i] <- xi_tmp$xi
      xi_acc_rate[i] <- xi_tmp$xi_acc
      prev_prop_xi_var <- xi_tmp$prop_xi_var

      # sigma2_S: holds omega fixed, so only the bGEV likelihood enters.
      sigma2_S_tmp <-
        sample_sigma2_S_base(
          sigma2_S[i - 1],
          omega,
          eta[i],
          nu[i],
          xi[i],
          y,
          stations = stations,
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
        sample_rho_S_base(
          rho_S[i - 1],
          omega,
          obs_coords,
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

      # Reconstruct S at the fully updated (sigma2_S[i], omega)
      S_k <- sqrt(sigma2_S[i]) * omega
      S_n[, i] <- S_k

      if (!(is.null(pred_coords))) {
        S_pred[, i] <-
          get_S_predlocs(rho_S[i], sigma2_S[i], S_n[, i], pred_coords, obs_coords, dist_11 = d11)
      }
    }

    set_priors <- list(
      eta = paste0("n(", PRIOR_ETA_MEAN, ",", PRIOR_ETA_VAR, ")"),
      nu = paste0("n(", PRIOR_NU_MEAN, ",", PRIOR_NU_VAR, ")"),
      xi = paste0("logn(", PRIOR_XI_MEAN, ",", PRIOR_XI_VAR, ")"),
      sigma2_S = paste0("IG(", PRIOR_SIGMA2_S_A, ",", PRIOR_SIGMA2_S_B, ")"),
      rho_S = paste0("Gamma(", PRIOR_RHO_S_A, ",", PRIOR_RHO_S_B, ")")
    )

    tmp_results <- list(
      eta = matrix(eta, nrow = 1),
      nu = matrix(nu, nrow = 1),
      xi = xi,
      S_n = S_n,
      S_pred = S_pred,
      sigma2_S = sigma2_S,
      rho_S = rho_S,
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
