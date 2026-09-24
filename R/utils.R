# source(here::here("R/constants.R"))

# Adaptive proposal scale
get_prop_var <- function(
    iter,
    prev_var,
    default_var,
    full_acc_rate,
    max_high,
    target = 0.44,
    window = 200) {
  if (is.null(iter)) {
    return(default_var)
  }
  adapt_until <- getOption("epsm.adapt_until", BURNIN)
  if (iter >= adapt_until) {
    return(prev_var)
  }
  if (is.null(full_acc_rate) || length(full_acc_rate) < 20) {
    return(prev_var)
  }
  recent <- mean(utils::tail(full_acc_rate, window), na.rm = TRUE)
  step <- min(0.5, 5 / sqrt(iter))
  prop_var <- prev_var * exp(step * (recent - target))
  return(min(max(prop_var, 1e-3), max_high))
}


make_prior_list <- function() {
  set_priors <- list(
    lambda = paste0("Gamma(", PRIOR_LAMBDA_A, ",", PRIOR_LAMBDA_B, ")"),
    eta = paste0("n(", PRIOR_ETA_MEAN, ",", PRIOR_ETA_VAR, ")"),
    nu = paste0("n(", PRIOR_NU_MEAN, ",", PRIOR_NU_VAR, ")"),
    xi = paste0("logn(", PRIOR_XI_MEAN, ",", PRIOR_XI_VAR, ")"),
    sigma2_S = paste0("IG(", PRIOR_SIGMA2_S_A, ",", PRIOR_SIGMA2_S_B, ")"),
    sigma2_T = paste0("IG(", PRIOR_SIGMA2_T_A, ",", PRIOR_SIGMA2_T_B, ")"),
    rho_S = paste0("Gamma(", PRIOR_RHO_S_A, ",", PRIOR_RHO_S_B, ")"),
    rho_T = paste0("Gamma(", PRIOR_RHO_T_A, ",", PRIOR_RHO_T_B, ")"),
    beta = paste0("n(", PRIOR_BETA_MEAN, ",", PRIOR_BETA_VAR, ")")
  )
  return(set_priors)
}

## Transform chains


get_chain_elements <- function(params, chain, burnin, thin = 1) {
  # Iterations to retain
  nsims <- length(chain$xi)
  rows_to_get <- seq(burnin, nsims, by = thin)

  merge_params <- function(p, chain) {
    if (
      p %in%
        c("lambda_star", "xi", "sigma2_W", "sigma2_S", "rho_W", "rho_S", "beta")
    ) {
      vec_pars <- matrix(
        chain[[p]][rows_to_get],
        ncol = 1
      )

      colnames(vec_pars) <- p
    } else if (p %in% c("eta", "nu", "S_n", "S_pred", "W_pred")) {
      vec_pars <- t(chain[[p]][, rows_to_get, drop = FALSE])

      colnames(vec_pars) <- paste0(
        p,
        ".",
        seq_len(ncol(vec_pars))
      )
    } else if (p == "W_n") {
      vec_pars <- purrr::map2_dfc(
        unname(chain$W_n),
        names(chain$W_n),
        ~ {
          # Select iterations FIRST, then transpose
          out <- t(.x[, rows_to_get, drop = FALSE])

          colnames(out) <- paste0(
            "W_n_",
            .y,
            ".",
            seq_len(ncol(out))
          )

          out
        }
      )

      vec_pars <- as.matrix(vec_pars)
    } else {
      stop("Unknown parameter: ", p)
    }

    vec_pars
  }

  # Build each parameter separately and column-bind them
  mcmc_elements <- do.call(
    cbind,
    lapply(params, merge_params, chain = chain)
  )

  mcmc_elements
}


#' Convert MCMC output to a coda mcmc.list
#'
#' @param res List of fitted chains (outputs of the `*_mcmc` functions).
#' @param burnin First iteration to keep.
#' @param thin Thinning interval.
#' @param params Parameter names to keep; default all sampled parameters.
#' @return A [coda::mcmc.list()].
#' @export
make_mcmc_list_obj <- function(res, burnin = 0, thin = 1, params = NULL) {
  # Note: res must be a list of chains.

  if (is.null(params)) {
    params <- names(res[[1]])[
      !(names(res[[1]])) %in%
        c("k", "pred_coords", "other_summaries", "sim_data")
    ]
  }
  all_elements <- lapply(
    res,
    function(mcmc, param) {
      coda::mcmc(get_chain_elements(param, mcmc, burnin, thin), start = 1)
    },
    params
  )
  return(coda::mcmc.list(all_elements))
}
