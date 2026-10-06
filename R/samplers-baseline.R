# source(here::here("R/bGEV/bGEVcode.R"))

#### Cholesky helpers

chol_field_factor <- function(coords, rho, jitter = 1e-6) {
  R <- exp_cor(fields::rdist(coords), range = rho) +
    diag(jitter, nrow(as.matrix(coords)))
  L <- FastGP::rcppeigen_get_chol(R)
  if (!all(is.finite(L))) {
    L <- t(chol(R))
  }
  return(L)
}

whiten_chol <- function(L, x) {
  return(forwardsolve(L, as.numeric(x)))
}

quad_chol <- function(L, x, y = NULL) {
  zx <- whiten_chol(L, x)
  if (is.null(y)) {
    return(sum(zx * zx))
  }
  return(sum(zx * whiten_chol(L, y)))
}

logdet_chol <- function(L) {
  return(2 * sum(log(diag(L))))
}

#### Joint elliptical slice sampler for (eta, omega)
lik_eta_omega_base <- function(eta, omega, sigma2_S, nu, xi, y, stations) {
  S_natural <- sqrt(sigma2_S) * omega
  q <- eta + S_natural[stations]
  sb <- rep(exp(nu), length(y))
  xi_vec <- rep(xi, length(y))
  sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
}

sample_eta_S_ess_base <- function(
    eta_cur,
    omega_cur,
    sigma2_S,
    rho_S,
    nu,
    xi,
    y,
    stations,
    coords,
    prior_eta_mean = 0,
    prior_eta_var = 10,
    n_ess_passes = 5,
    chol_R = NULL) {
  k <- length(omega_cur)
  if (is.null(chol_R)) {
    chol_R <- chol_field_factor(coords, rho_S)
  }

  sd_eta <- sqrt(prior_eta_var)
  state_cur <- c(eta_cur - prior_eta_mean, omega_cur)

  lik_of <- function(state) {
    eta_val <- state[1] + prior_eta_mean
    omega_val <- state[-1]
    lik_eta_omega_base(eta_val, omega_val, sigma2_S, nu, xi, y, stations)
  }

  lik_cur <- lik_of(state_cur)

  for (pass in seq_len(n_ess_passes)) {
    eps_eta <- rnorm(1, 0, sd_eta)
    eps_omega <- as.vector(chol_R %*% rnorm(k))
    ellipse_ess <- c(eps_eta, eps_omega)

    logy <- lik_cur + log(runif(1))

    theta <- runif(1, 0, 2 * pi)
    theta_min <- theta - 2 * pi
    theta_max <- theta

    repeat {
      state_prop <- state_cur * cos(theta) + ellipse_ess * sin(theta)
      lik_prop <- lik_of(state_prop)
      if (is.finite(lik_prop) && (lik_prop > logy)) {
        state_cur <- state_prop
        lik_cur <- lik_prop
        break
      }
      if (theta > 0) {
        theta_max <- theta
      } else {
        theta_min <- theta
      }
      if ((theta_max - theta_min) < 1e-12) {
        break
      }
      theta <- runif(1, theta_min, theta_max)
    }
  }

  return(list(eta = state_cur[1] + prior_eta_mean, omega = state_cur[-1]))
}

#### nu
sample_nu_base <- function(
    nu_cur,
    eta,
    xi,
    S_n,
    y,
    prior_nu_mean = 0,
    prior_nu_var = 1,
    iter = NULL,
    prev_prop_nu_var = NULL,
    full_acc_rate = NULL) {
  nu_acc <- 0
  prop_nu_var <- get_prop_var(
    iter, prev_prop_nu_var,
    default_var = 0.5, full_acc_rate, max_high = 3
  )
  nu_prop <- rnorm(1, nu_cur, prop_nu_var)
  log_mh_ratio <-
    target_nu_base(nu_prop, eta, xi, S_n, y, prior_nu_mean, prior_nu_var) -
    target_nu_base(nu_cur, eta, xi, S_n, y, prior_nu_mean, prior_nu_var)

  if (runif(1) < exp(log_mh_ratio)) {
    nu_acc <- 1
    nu_cur <- nu_prop
  }
  return(list(nu = nu_cur, nu_acc = nu_acc, prop_nu_var = prop_nu_var))
}

target_nu_base <- function(nu_tmp, eta, xi, S_n, y, prior_nu_mean, prior_nu_var) {
  q <- eta + S_n
  sb <- rep(exp(nu_tmp), length(y))
  xi_vec <- rep(xi, length(y))
  t1 <- sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
  t2 <- (-1 / 2) * (nu_tmp - prior_nu_mean)^2 / prior_nu_var
  return(t1 + t2)
}

#### bGEV-only likelihood
lik_y_bgev_base <- function(S_vals, eta, nu, xi, y, stations) {
  S_vals <- as.numeric(S_vals)
  q <- eta + S_vals[stations]
  sb <- rep(exp(nu), length(y))
  sum(dbgev2_vec(y, q, sb, rep(xi, length(y)), log = TRUE))
}

#### sigma2_S
sample_sigma2_S_base <- function(
    sigma2_S_cur,
    omega,
    eta,
    nu,
    xi,
    y,
    stations,
    sigma2_S_prior_a,
    sigma2_S_prior_b,
    iter = NULL,
    prev_prop_sigma2_S_var = NULL,
    full_acc_rate = NULL) {
  sigma2_S_acc <- 0
  omega <- as.numeric(omega)

  prop_sigma2_S_var <- get_prop_var(
    iter, prev_prop_sigma2_S_var,
    default_var = 0.25, full_acc_rate, max_high = 1.0
  )
  sigma2_S_prop <- exp(rnorm(1, log(sigma2_S_cur), prop_sigma2_S_var))

  S_cur <- sqrt(sigma2_S_cur) * omega
  S_prop <- sqrt(sigma2_S_prop) * omega

  log_mh_ratio <-
    lik_y_bgev_base(S_prop, eta, nu, xi, y, stations) -
    lik_y_bgev_base(S_cur, eta, nu, xi, y, stations) -
    sigma2_S_prior_a * (log(sigma2_S_prop) - log(sigma2_S_cur)) -
    sigma2_S_prior_b * (1 / sigma2_S_prop - 1 / sigma2_S_cur)

  if (is.finite(log_mh_ratio) && (log(runif(1)) < log_mh_ratio)) {
    sigma2_S_acc <- 1
    sigma2_S_cur <- sigma2_S_prop
  }
  return(list(
    sigma2_S = sigma2_S_cur,
    sigma2_S_acc = sigma2_S_acc,
    prop_sigma2_S_var = prop_sigma2_S_var
  ))
}

#### rho_S
rho_mh_logratio_base <- function(
    rho_cur, rho_prop, gp_vec, rho_prior_a, rho_prior_b, coords_cur, chol_R_cur = NULL) {
  U_prop <- chol_field_factor(coords_cur, rho_prop)
  if (is.null(chol_R_cur)) {
    chol_R_cur <- chol_field_factor(coords_cur, rho_cur)
  }
  term1 <- (-1 / 2) * (logdet_chol(U_prop) - logdet_chol(chol_R_cur))
  term2 <- rho_prior_a * (log(rho_prop) - log(rho_cur))
  term3 <-
    -(1 / 2) * (quad_chol(U_prop, gp_vec) - quad_chol(chol_R_cur, gp_vec)) -
    rho_prior_b * (rho_prop - rho_cur)
  return(term1 + term2 + term3)
}

sample_rho_S_base <- function(
    rho_cur,
    omega,
    obs_coords,
    rho_prior_a = 1000,
    rho_prior_b = 5000,
    iter = NULL,
    prev_prop_rho_var = 0.35,
    full_acc_rate = NULL,
    chol_R_cur = NULL) {
  rho_acc <- 0
  prop_rho_var <- get_prop_var(
    iter, prev_prop_rho_var,
    default_var = 0.35, full_acc_rate, max_high = 1.0
  )
  rho_prop <- exp(rnorm(1, log(rho_cur), prop_rho_var))

  log_mh_ratio_rho <- rho_mh_logratio_base(
    rho_cur, rho_prop, omega, rho_prior_a, rho_prior_b, obs_coords, chol_R_cur
  )

  if (is.finite(log_mh_ratio_rho) && (log(runif(1)) < log_mh_ratio_rho)) {
    rho_acc <- 1
    rho_cur <- rho_prop
  }
  return(list(rho_S = rho_cur, rho_S_acc = rho_acc, prop_rho_var = prop_rho_var))
}

#### xi
sample_xi_base <- function(
    xi_cur,
    eta,
    nu,
    S_n,
    y,
    stations,
    prior_xi_mean = -1,
    prior_xi_var = 1,
    iter = NULL,
    prev_prop_xi_var = NULL,
    full_acc_rate = NULL) {
  xi_acc <- 0
  prop_xi_var <- 0.7
  xi_prop <- rlnorm(1, log(xi_cur), prop_xi_var)
  num <- target_xi_base(xi_prop, eta, nu, S_n, y, stations, prior_xi_mean, prior_xi_var)
  den <- target_xi_base(xi_cur, eta, nu, S_n, y, stations, prior_xi_mean, prior_xi_var)
  log_mh_ratio <- num - den - log(xi_cur) + log(xi_prop)
  if (is.na(exp(log_mh_ratio))) {
    return(list(xi = xi_cur, xi_acc = xi_acc, prop_xi_var = prop_xi_var))
  }
  if (runif(1) < exp(log_mh_ratio)) {
    xi_acc <- 1
    xi_cur <- xi_prop
  }
  return(list(xi = xi_cur, xi_acc = xi_acc, prop_xi_var = prop_xi_var))
}

target_xi_base <- function(xi_tmp, eta, nu, S_n, y, stations, prior_xi_mean, prior_xi_var) {
  q <- eta + S_n[stations]
  sb <- rep(exp(nu), length(y))
  xi_vec <- rep(xi_tmp, length(y))
  t1 <- sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
  t2 <- dlnorm(xi_tmp, prior_xi_mean, sdlog = sqrt(prior_xi_var), log = TRUE)
  return(t1 + t2)
}

# ===========================================================================
# Per-year (baseline-timevar.R) versions of sigma2_S, rho_S, and xi
# ===========================================================================

# q_j(s_i) = eta_j + S_j(s_i), assembled per-observation. flat_idx[i] is the
# position of observation i's (year, local-station) pair within
# unlist(S_n_natural_list)
make_q_vector_base <- function(idx_years, eta, S_n_natural_list, flat_idx) {
  eta[idx_years] + unlist(S_n_natural_list, use.names = FALSE)[flat_idx]
}

lik_y_bgev_base_timevar <- function(S_vals_list, eta, nu, xi, y, flat_idx, idx_years) {
  q <- make_q_vector_base(idx_years, eta, S_vals_list, flat_idx)
  sb <- exp(nu[idx_years])
  sum(dbgev2_vec(y, q, sb, rep(xi, length(y)), log = TRUE))
}

sample_sigma2_S_base_timevar <- function(
    sigma2_S_cur,
    omega_list,
    eta,
    nu,
    xi,
    y,
    flat_idx,
    idx_years,
    sigma2_S_prior_a,
    sigma2_S_prior_b,
    iter = NULL,
    prev_prop_sigma2_S_var = NULL,
    full_acc_rate = NULL) {
  sigma2_S_acc <- 0

  prop_sigma2_S_var <- get_prop_var(
    iter, prev_prop_sigma2_S_var,
    default_var = 0.25, full_acc_rate, max_high = 1.0
  )
  sigma2_S_prop <- exp(rnorm(1, log(sigma2_S_cur), prop_sigma2_S_var))

  S_cur_list <- lapply(omega_list, function(x) sqrt(sigma2_S_cur) * x)
  S_prop_list <- lapply(omega_list, function(x) sqrt(sigma2_S_prop) * x)

  log_mh_ratio <-
    lik_y_bgev_base_timevar(S_prop_list, eta, nu, xi, y, flat_idx, idx_years) -
    lik_y_bgev_base_timevar(S_cur_list, eta, nu, xi, y, flat_idx, idx_years) -
    sigma2_S_prior_a * (log(sigma2_S_prop) - log(sigma2_S_cur)) -
    sigma2_S_prior_b * (1 / sigma2_S_prop - 1 / sigma2_S_cur)

  if (is.finite(log_mh_ratio) && (log(runif(1)) < log_mh_ratio)) {
    sigma2_S_acc <- 1
    sigma2_S_cur <- sigma2_S_prop
  }
  return(list(
    sigma2_S = sigma2_S_cur,
    sigma2_S_acc = sigma2_S_acc,
    prop_sigma2_S_var = prop_sigma2_S_var
  ))
}

# omega_cur_list
rho_mh_logratio_base_timevar <- function(
    rho_cur, R_cur, R_cur_inv, rho_prop, omega_cur_list,
    rho_prior_a, rho_prior_b, coords_cur, stations_yearly_index) {
  n_years <- length(omega_cur_list)

  logdet_prop <- numeric(n_years)
  logdet_cur <- numeric(n_years)
  Q_prop <- 0
  Q_cur <- 0

  for (yr in seq_len(n_years)) {
    dist_yr <- fields::rdist(coords_cur[stations_yearly_index[[yr]], , drop = FALSE])
    R_prop_yr <- exp_cor(dist_yr, range = rho_prop) + diag(1e-8, nrow(dist_yr))
    L_prop_yr <- FastGP::rcppeigen_get_chol(R_prop_yr)
    R_prop_inv_yr <- chol2inv(t(L_prop_yr))

    L_cur_yr <- FastGP::rcppeigen_get_chol(R_cur[[yr]] + diag(1e-8, nrow(R_cur[[yr]])))

    logdet_prop[yr] <- 2 * sum(log(diag(L_prop_yr)))
    logdet_cur[yr] <- 2 * sum(log(diag(L_cur_yr)))

    Q_prop <- Q_prop + t(omega_cur_list[[yr]]) %*% R_prop_inv_yr %*% omega_cur_list[[yr]]
    Q_cur <- Q_cur + t(omega_cur_list[[yr]]) %*% R_cur_inv[[yr]] %*% omega_cur_list[[yr]]
  }

  term1 <- (-1 / 2) * sum(logdet_prop - logdet_cur)
  term2 <- rho_prior_a * (log(rho_prop) - log(rho_cur))
  term3 <- -(1 / 2) * as.numeric(Q_prop - Q_cur) - rho_prior_b * (rho_prop - rho_cur)

  return(term1 + term2 + term3)
}

sample_rho_S_base_timevar <- function(
    rho_cur, R_cur, R_cur_inv, omega_cur_list, obs_coords, stations_yearly_index,
    rho_prior_a, rho_prior_b) {
  rho_acc <- 0
  rho_prop <- stats::rlnorm(1, log(rho_cur), 0.1)

  log_mh_ratio_rho <- rho_mh_logratio_base_timevar(
    rho_cur, R_cur, R_cur_inv, rho_prop, omega_cur_list,
    rho_prior_a, rho_prior_b, obs_coords, stations_yearly_index
  )

  if (runif(1) < exp(log_mh_ratio_rho)) {
    rho_acc <- 1
    rho_cur <- ifelse(rho_prop < 0.001, rho_cur, rho_prop)
  }
  return(list(rho_S = rho_cur, rho_S_acc = rho_acc))
}

sample_xi_base_timevar <- function(
    xi_cur,
    eta,
    nu,
    S_n_list,
    y,
    flat_idx,
    idx_years,
    prior_xi_mean = -1,
    prior_xi_var = 1,
    iter = NULL,
    prev_prop_xi_var = NULL,
    full_acc_rate = NULL) {
  xi_acc <- 0
  prop_xi_var <- 0.7
  xi_prop <- rlnorm(1, log(xi_cur), prop_xi_var)
  num <- target_xi_base_timevar(
    xi_prop, eta, nu, S_n_list, y, flat_idx, idx_years, prior_xi_mean, prior_xi_var
  )
  den <- target_xi_base_timevar(
    xi_cur, eta, nu, S_n_list, y, flat_idx, idx_years, prior_xi_mean, prior_xi_var
  )
  log_mh_ratio <- num - den - log(xi_cur) + log(xi_prop)
  if (is.na(exp(log_mh_ratio))) {
    return(list(xi = xi_cur, xi_acc = xi_acc, prop_xi_var = prop_xi_var))
  }
  if (runif(1) < exp(log_mh_ratio)) {
    xi_acc <- 1
    xi_cur <- xi_prop
  }
  return(list(xi = xi_cur, xi_acc = xi_acc, prop_xi_var = prop_xi_var))
}

target_xi_base_timevar <- function(
    xi_tmp, eta, nu, S_n_list, y, flat_idx, idx_years, prior_xi_mean, prior_xi_var) {
  q <- make_q_vector_base(idx_years, eta, S_n_list, flat_idx)
  sb <- exp(nu[idx_years])
  xi_vec <- rep(xi_tmp, length(y))
  t1 <- sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
  t2 <- dlnorm(xi_tmp, prior_xi_mean, sdlog = sqrt(prior_xi_var), log = TRUE)
  return(t1 + t2)
}
