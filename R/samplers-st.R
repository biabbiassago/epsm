# source(here::here("R/bGEV/bGEVcode.R"))

# q(s_i) = eta_j + W_j(s_i) + S(s_i), assembled per-observation,
make_q_vector <- function(idx_years, eta, S_n, W_n, stations, flat_idx) {
  eta[idx_years] + unlist(W_n, use.names = FALSE)[flat_idx] + S_n[stations]
}
oneyear_likelihood_bGEV_timevar <- function(y, eta, S_n, W_n, nu, xi) {
  q <- eta + S_n + W_n
  sb <- rep(exp(nu), length(y))
  xi <- rep(xi, length(y))

  tmp <-
    cbind(
      "y" = y,
      "q" = q,
      "sb" = sb,
      "xi" = xi
    )
  t1 <- sum(dbgev2_vec(
    tmp[, "y"],
    tmp[, "q"],
    tmp[, "sb"],
    tmp[, "xi"],
    log = TRUE
  ))
  return(t1)
}

#### Poisson Process Part ####
sample_lambda_star <-
  function(k_cur, lambda_prior_a, lambda_prior_b, area_B = 1) {
    lambda_star <- rgamma(1, lambda_prior_a + k_cur, lambda_prior_b + area_B)
    return(lambda_star)
  }

sample_all_coords <-
  function(
      lambda_star,
      beta,
      rho_S,
      omega_k,
      obs_coords,
      all_coords_prev,
      owin = owin(),
      area_B = 1) {
    # Simulate a homogeneous pp with intensity lambda star
    k <- dim(all_coords_prev)[1]
    k_star <- rpois(1, lambda_star * area_B)
    if (k_star < 1) {
      return(list("all_coords" = obs_coords, "omega_x_tilde" = numeric(0)))
    }
    pp <- spatstat.random::runifpoint(k_star, owin = owin)
    coords_pp <- spatstat.geom::coords(pp)

    big_coords <- rbind(all_coords_prev, coords_pp)
    R_full <- exp_cor(fields::rdist(big_coords), range = rho_S) +
      diag(1e-6, k + k_star)
    L_full <- FastGP::rcppeigen_get_chol(R_full)

    i1 <- seq_len(k)
    i2 <- k + seq_len(k_star)
    omega_1 <- forwardsolve(L_full[i1, i1, drop = FALSE], as.numeric(omega_k))
    omega_x_tilde_all <- as.numeric(
      L_full[i2, i1, drop = FALSE] %*%
        omega_1 +
        L_full[i2, i2, drop = FALSE] %*% rnorm(k_star)
    )

    probs <- pnorm(-beta * omega_x_tilde_all)
    idx_keep <- rbinom(k_star, 1, probs) == 1

    disc_coords <- coords_pp[idx_keep, ]
    omega_x_tilde <- omega_x_tilde_all[idx_keep]

    all_coords <- rbind(obs_coords, disc_coords)
    return(list("all_coords" = all_coords, "omega_x_tilde" = omega_x_tilde))
  }

#### Cholesky helpers ####
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

### Spatial field + eta: per-year joint elliptical slice sampler (W held fixed) ####

sample_eta_j_S_ess_timevar <- function(
    eta_full,
    j,
    omega_k_cur,
    beta,
    sigma2_S,
    rho_S,
    W_n,
    nu,
    xi,
    y,
    littlen,
    all_coords_new,
    stations,
    flat_idx,
    idx_years,
    prior_eta_mean = 0,
    prior_eta_var = 10,
    n_ess_passes = 3,
    chol_R = NULL) {
  k_new <- dim(all_coords_new)[1]
  if (length(omega_k_cur) != k_new) {
    stop(
      "sample_eta_j_S_ess_timevar: field length does not match the number of locations"
    )
  }

  if (is.null(chol_R)) {
    chol_R <- chol_field_factor(all_coords_new, rho_S)
  }

  sd_eta <- sqrt(prior_eta_var)
  eta_j_cur <- eta_full[j] - prior_eta_mean
  state_cur <- c(eta_j_cur, omega_k_cur)

  lik_of <- function(state) {
    eta_vec <- eta_full
    eta_vec[j] <- state[1] + prior_eta_mean
    omega_vec <- state[-1]
    lik_eta_omega_timevar(
      eta_vec,
      omega_vec,
      W_n,
      beta,
      sigma2_S,
      nu,
      xi,
      y,
      littlen,
      k_new,
      stations,
      flat_idx,
      idx_years
    )
  }

  lik_cur <- lik_of(state_cur)

  for (pass in seq_len(n_ess_passes)) {
    eps_eta <- rnorm(1, 0, sd_eta)
    eps_omega <- as.vector(chol_R %*% rnorm(k_new))
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

  eta_j_new <- state_cur[1] + prior_eta_mean
  omega_k_new <- state_cur[-1]
  return(list(eta_j = eta_j_new, omega_k = omega_k_new))
}

# Joint log-target for (eta, omega), W held fixed.
lik_eta_omega_timevar <-
  function(
      eta_vec,
      omega_k,
      W_n,
      beta,
      sigma2_S,
      nu,
      xi,
      y,
      littlen,
      k,
      stations,
      flat_idx,
      idx_years) {
    I_k <- c(rep(1, littlen), rep(-1, k - littlen))
    t1 <- sum(pnorm(beta * I_k * omega_k, log = TRUE))

    sigma_S <- sqrt(sigma2_S)
    S_n_natural <- sigma_S * omega_k
    q <- make_q_vector(
      idx_years,
      eta_vec,
      S_n_natural,
      W_n,
      stations,
      flat_idx
    )
    sb <- exp(nu[idx_years])
    xi_vec <- rep(xi, length(y))

    t2 <- sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
    return(t1 + t2)
  }

# bGEV part of the likelihood only (no point-process term), W held fixed and
# added back in via make_q_vector. Used by sample_sigma2_S, which never
# touches the point process under partial whitening.
lik_y_bgev_timevar <- function(
    S_vals,
    eta,
    W_n,
    nu,
    xi,
    y,
    stations,
    flat_idx,
    idx_years) {
  S_vals <- as.numeric(S_vals)
  q <- make_q_vector(idx_years, eta, S_vals, W_n, stations, flat_idx)
  sb <- exp(nu[idx_years])
  sum(dbgev2_vec(y, q, sb, rep(xi, length(y)), log = TRUE))
}

#### bGEV Part

sample_nu_timevar <-
  function(
      nu_cur,
      eta,
      xi,
      S_n,
      W_n,
      y,
      stations = NULL,
      varying_range,
      prior_nu_mean = 0,
      prior_nu_var = 1,
      iter = NULL,
      prev_prop_nu_var = NULL,
      full_acc_rate = NULL) {
    nu_acc <- 0
    prop_nu_var <- get_prop_var(
      iter,
      prev_prop_nu_var,
      default_var = 0.5,
      full_acc_rate,
      max_high = 3
    )

    nu_prop <- rnorm(1, nu_cur, prop_nu_var)
    log_mh_ratio <-
      target_nu_timevar(
        nu_prop,
        eta,
        xi,
        S_n,
        W_n,
        y,
        stations,
        varying_range,
        prior_nu_mean,
        prior_nu_var
      ) -
      target_nu_timevar(
        nu_cur,
        eta,
        xi,
        S_n,
        W_n,
        y,
        stations,
        varying_range,
        prior_nu_mean,
        prior_nu_var
      )

    if (runif(1) < exp(log_mh_ratio)) {
      nu_acc <- 1
      nu_cur <- nu_prop
    }
    return(list(nu = nu_cur, nu_acc = nu_acc, prop_nu_var = prop_nu_var))
  }
target_nu_timevar <-
  function(
      nu_tmp,
      eta,
      xi,
      S_n,
      W_n,
      y,
      stations = NULL,
      varying_range,
      prior_nu_mean,
      prior_nu_var) {
    t1 <- oneyear_likelihood_bGEV_timevar(y, eta, S_n, W_n, nu_tmp, xi)
    t2 <- (-1 / 2) * (nu_tmp - prior_nu_mean)^2 / prior_nu_var
    return(t1 + t2)
  }

### W(s): whitened elliptical slice sampler ####

sample_omega_W_ess_timevar <- function(
    omega_W_cur,
    R_W,
    sigma2_W,
    eta,
    nu,
    xi,
    y,
    S_n) {
  epsilon <- 1e-6 # Small jitter for numerical stability
  SIGMA <- R_W + diag(epsilon, length(omega_W_cur)) # unit sill: no sigma2_W here
  ellipse_ess <- as.vector(FastGP::rcpp_rmvnorm(
    1,
    S = SIGMA,
    mu = rep(0, length(omega_W_cur))
  ))
  u <- runif(1)

  sigma_W <- sqrt(sigma2_W)
  lik_of <- function(omega_W) {
    oneyear_likelihood_bGEV_timevar(y, eta, S_n, sigma_W * omega_W, nu, xi)
  }

  logy <- lik_of(omega_W_cur) + log(u)

  theta <- runif(1, 0, 2 * pi)
  theta_min <- theta - 2 * pi
  theta_max <- theta

  ess_omega_w_targ <- function(theta) {
    omega_W_new <- omega_W_cur * cos(theta) + ellipse_ess * sin(theta)
    tmp_targ <- lik_of(omega_W_new)
    return(list("s" = omega_W_new, "t" = tmp_targ))
  }
  ess_step <- ess_omega_w_targ(theta)
  tmp_targ <- ess_step$t

  while_iter <- 1
  while (tmp_targ < logy) {
    while_iter <- while_iter + 1
    if (while_iter > 10000) {
      stop("Elliptical Slice Sampler for omega_W not converging")
    }

    if (theta < 0) {
      theta_min <- theta
    } else {
      theta_max <- theta
    }

    theta <- runif(1, theta_min, theta_max)
    ess_step <- ess_omega_w_targ(theta)
    tmp_targ <- ess_step$t
  }
  return(ess_step$s)
}

# sigma2_W under whitening
# while holding W_n fixed
sample_sigma2_W_timevar <- function(
    sigma2_W_cur,
    omega_W_list,
    eta,
    S_n,
    nu,
    xi,
    y,
    stations,
    flat_idx,
    idx_years,
    sigma2_W_prior_a,
    sigma2_W_prior_b,
    prop_sd = 0.25) {
  sigma2_W_acc <- 0
  sigma2_W_prop <- exp(rnorm(1, log(sigma2_W_cur), prop_sd))

  W_n_cur <- lapply(omega_W_list, function(x) sqrt(sigma2_W_cur) * x)
  W_n_prop <- lapply(omega_W_list, function(x) sqrt(sigma2_W_prop) * x)

  log_mh_ratio <-
    lik_y_bgev_timevar(
      S_n,
      eta,
      W_n_prop,
      nu,
      xi,
      y,
      stations,
      flat_idx,
      idx_years
    ) -
    lik_y_bgev_timevar(
      S_n,
      eta,
      W_n_cur,
      nu,
      xi,
      y,
      stations,
      flat_idx,
      idx_years
    ) -
    sigma2_W_prior_a * (log(sigma2_W_prop) - log(sigma2_W_cur)) -
    sigma2_W_prior_b * (1 / sigma2_W_prop - 1 / sigma2_W_cur)

  if (is.finite(log_mh_ratio) && (log(runif(1)) < log_mh_ratio)) {
    sigma2_W_acc <- 1
    sigma2_W_cur <- sigma2_W_prop
  }
  return(list(sigma2_W = sigma2_W_cur, sigma2_W_acc = sigma2_W_acc))
}

sample_xi_timevar <-
  function(
      xi_cur,
      eta,
      nu,
      S_n,
      W_n,
      y,
      stations,
      flat_idx,
      idx_years,
      prior_xi_mean = -1,
      prior_xi_var = 1,
      varying_range,
      iter = NULL,
      prev_prop_xi_var = NULL,
      full_acc_rate = NULL) {
    xi_acc <- 0
    prop_xi_var <- 0.7
    xi_prop <- rlnorm(1, log(xi_cur), prop_xi_var)
    num <-
      target_xi_timevar(
        xi_prop,
        eta,
        nu,
        S_n,
        W_n,
        y,
        stations,
        flat_idx,
        idx_years,
        prior_xi_mean,
        prior_xi_var,
        varying_range = varying_range
      )
    den <-
      target_xi_timevar(
        xi_cur,
        eta,
        nu,
        S_n,
        W_n,
        y,
        stations,
        flat_idx,
        idx_years,
        prior_xi_mean,
        prior_xi_var,
        varying_range = varying_range
      )
    log_mh_ratio <-
      num - den - log(xi_cur) + log(xi_prop)
    if (is.na(exp(log_mh_ratio))) {
      return(list(xi = xi_cur, xi_acc = xi_acc))
    }
    if (runif(1) < exp(log_mh_ratio)) {
      xi_acc <- 1
      xi_cur <- xi_prop
    }
    return(list(xi = xi_cur, xi_acc = xi_acc, prop_xi_var = prop_xi_var))
  }
target_xi_timevar <-
  function(
      xi_tmp,
      eta,
      nu,
      S_n,
      W_n,
      y,
      stations,
      flat_idx,
      idx_years,
      prior_xi_mean,
      prior_xi_var,
      varying_range) {
    # Remember xi must always be positive.
    q <- make_q_vector(idx_years, eta, S_n, W_n, stations, flat_idx)

    if (varying_range == T) {
      sb <- exp(nu[idx_years] + T_n[stations])
      stop("Sorry currently not handled")
    } else if (varying_range == F) {
      sb <- exp(nu[idx_years])
    }
    xi <- rep(xi_tmp, length(y))

    tmp <-
      cbind(
        "y" = y,
        "q" = q,
        "sb" = sb,
        "xi" = xi
      )

    t1 <- sum(dbgev2_vec(
      tmp[, "y"],
      tmp[, "q"],
      tmp[, "sb"],
      tmp[, "xi"],
      log = TRUE
    ))
    t2 <-
      dlnorm(xi_tmp, prior_xi_mean, sdlog = sqrt(prior_xi_var), log = T)
    return(t1 + t2)
  }

#### Spatial Hyper-parameters Part (shared field S) ####

# sigma2_S under partial whitening
sample_sigma2_S_timevar <-
  function(sigma2_S_cur,
           omega_k,
           eta,
           W_n,
           nu,
           xi,
           y,
           stations,
           flat_idx,
           idx_years,
           sigma2_S_prior_a,
           sigma2_S_prior_b,
           iter = NULL,
           prev_prop_sigma2_S_var = NULL,
           full_acc_rate = NULL) {
    sigma2_S_acc <- 0
    omega_k <- as.numeric(omega_k)

    prop_sigma2_S_var <- get_prop_var(
      iter,
      prev_prop_sigma2_S_var,
      default_var = 0.25,
      full_acc_rate,
      max_high = 1.0
    )

    sigma2_S_prop <- exp(rnorm(1, log(sigma2_S_cur), prop_sigma2_S_var))

    S_cur <- sqrt(sigma2_S_cur) * omega_k
    S_prop <- sqrt(sigma2_S_prop) * omega_k

    log_mh_ratio <-
      lik_y_bgev_timevar(
        S_prop,
        eta,
        W_n,
        nu,
        xi,
        y,
        stations,
        flat_idx,
        idx_years
      ) -
      lik_y_bgev_timevar(
        S_cur,
        eta,
        W_n,
        nu,
        xi,
        y,
        stations,
        flat_idx,
        idx_years
      ) -
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

# rho_S
sample_rho_S <-
  function(
      rho_cur,
      omega_k,
      all_coords,
      rho_prior_a = 1000,
      rho_prior_b = 5000,
      iter = NULL,
      prev_prop_rho_var = 0.35,
      full_acc_rate = NULL,
      chol_R_cur = NULL) {
    rho_acc <- 0
    prop_rho_var <-
      get_prop_var(
        iter,
        prev_prop_rho_var,
        default_var = 0.35,
        full_acc_rate,
        max_high = 1.0
      )
    log_rho_prop <- stats::rnorm(1, mean = log(rho_cur), sd = prop_rho_var)
    rho_prop <- exp(log_rho_prop)

    log_mh_ratio_rho <-
      rho_mh_logratio(
        rho_cur = rho_cur,
        rho_prop = rho_prop,
        gp_vec = omega_k,
        chol_R_cur = chol_R_cur,
        rho_prior_a = rho_prior_a,
        rho_prior_b = rho_prior_b,
        coords_cur = all_coords
      )

    if (is.finite(log_mh_ratio_rho) && (log(runif(1)) < log_mh_ratio_rho)) {
      rho_acc <- 1
      rho_cur <- rho_prop
    }
    return(list(
      rho_S = rho_cur,
      rho_S_acc = rho_acc,
      prop_rho_var = prop_rho_var
    ))
  }

rho_mh_logratio <-
  function(
      rho_cur,
      rho_prop,
      gp_vec,
      rho_prior_a,
      rho_prior_b,
      coords_cur,
      chol_R_cur = NULL) {
    U_prop <- chol_field_factor(coords_cur, rho_prop)
    if (is.null(chol_R_cur)) {
      chol_R_cur <- chol_field_factor(coords_cur, rho_cur)
    }

    term1 <-
      (-1 / 2) * (logdet_chol(U_prop) - logdet_chol(chol_R_cur))

    term2 <- (rho_prior_a) * (log(rho_prop) - log(rho_cur))

    term3 <-
      -(1 / 2) *
      (quad_chol(U_prop, gp_vec) - quad_chol(chol_R_cur, gp_vec)) -
      (rho_prior_b) * (rho_prop - rho_cur)

    return(term1 + term2 + term3)
  }

#### Preferential Parameter Part ####

# beta
sample_beta <-
  function(
      beta_cur,
      omega_k,
      littlen,
      k,
      beta_prior_mean = 0,
      beta_prior_var = 1,
      iter = NULL,
      prev_prop_beta_var = NULL,
      full_acc_rate = NULL) {
    beta_acc <- 0
    prop_beta_var <- get_prop_var(
      iter,
      prev_prop_beta_var,
      default_var = 0.7,
      full_acc_rate,
      max_high = 2.0
    )
    beta_prop <- rnorm(1, beta_cur, prop_beta_var)

    I_k <- c(rep(1, littlen), rep(-1, k - littlen))
    log_mh_ratio <-
      target_beta_log(
        beta_prop,
        I_k,
        omega_k,
        beta_prior_mean,
        beta_prior_var
      ) -
      target_beta_log(beta_cur, I_k, omega_k, beta_prior_mean, beta_prior_var)

    if (runif(1) < exp(log_mh_ratio)) {
      beta_acc <- 1
      beta_cur <- beta_prop
    }
    return(list(
      beta = beta_cur,
      beta_acc = beta_acc,
      prop_beta_var = prop_beta_var
    ))
  }

target_beta_log <-
  function(beta, I_k, omega_k, beta_prior_mean, beta_prior_var) {
    return(
      sum(pnorm(beta * I_k * omega_k, log = TRUE)) +
        dnorm(beta, beta_prior_mean, sqrt(beta_prior_var), log = TRUE)
    )
  }

#### Year-specific random effect W(s):

# Sample rho_W
sample_rho_base_timevar <-
  function(
      rho_cur,
      R_W_cur,
      R_W_cur_inv,
      omega_W_cur,
      obs_coords,
      stations_yearly_index,
      rho_prior_a,
      rho_prior_b) {
    rho_acc <- 0
    rho_prop <- stats::rlnorm(1, log(rho_cur), 0.1)

    log_mh_ratio_rho <- rho_mh_logratio_timevar(
      rho_cur,
      R_W_cur,
      R_W_cur_inv,
      rho_prop,
      omega_W_cur,
      rho_prior_a,
      rho_prior_b,
      obs_coords,
      stations_yearly_index
    )

    if (runif(1) < exp(log_mh_ratio_rho)) {
      rho_acc <- 1
      rho_cur <- ifelse(rho_prop < 0.001, rho_cur, rho_prop)
    }
    return(list(rho_W = rho_cur, rho_W_acc = rho_acc))
  }

rho_mh_logratio_timevar <- function(
    rho_cur,
    R_cur,
    R_cur_inv,
    rho_prop,
    omega_W_cur,
    rho_prior_a,
    rho_prior_b,
    coords_cur,
    stations_yearly_index) {
  n_years <- length(omega_W_cur)

  logdet_prop <- numeric(n_years)
  logdet_cur <- numeric(n_years)
  Q_prop <- 0
  Q_cur <- 0

  for (yr in seq_len(n_years)) {
    dist_yr <- fields::rdist(coords_cur[
      stations_yearly_index[[yr]], ,
      drop = FALSE
    ])
    R_prop_yr <- exp_cor(dist_yr, range = rho_prop) + diag(1e-8, nrow(dist_yr))
    L_prop_yr <- FastGP::rcppeigen_get_chol(R_prop_yr)
    R_prop_inv_yr <- chol2inv(t(L_prop_yr))

    L_cur_yr <- FastGP::rcppeigen_get_chol(
      R_cur[[yr]] + diag(1e-8, nrow(R_cur[[yr]]))
    )

    logdet_prop[yr] <- 2 * sum(log(diag(L_prop_yr)))
    logdet_cur[yr] <- 2 * sum(log(diag(L_cur_yr)))

    Q_prop <- Q_prop +
      t(omega_W_cur[[yr]]) %*% R_prop_inv_yr %*% omega_W_cur[[yr]]
    Q_cur <- Q_cur +
      t(omega_W_cur[[yr]]) %*% R_cur_inv[[yr]] %*% omega_W_cur[[yr]]
  }

  term1 <- (-1 / 2) * sum(logdet_prop - logdet_cur)
  term2 <- (rho_prior_a) * (log(rho_prop) - log(rho_cur))
  term3 <- -(1 / 2) *
    as.numeric(Q_prop - Q_cur) -
    rho_prior_b * (rho_prop - rho_cur)

  return(term1 + term2 + term3)
}
