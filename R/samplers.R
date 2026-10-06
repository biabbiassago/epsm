# source(here::here("R/bGEV/bGEVcode.R"))

#### Poisson Process Part ####
sample_lambda_star <-
  function(k_cur, lambda_prior_a, lambda_prior_b, area_B = 1) {
    lambda_star <- rgamma(1, lambda_prior_a + k_cur, lambda_prior_b + area_B)
    return(lambda_star)
  }

sample_all_coords <-
  function(lambda_star,
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
    # gotta generate them in the actual area of interest.
    pp <- spatstat.random::runifpoint(k_star, win = owin)
    coords_pp <- spatstat.geom::coords(pp)

    # Draw the field at the new points conditional on the current field using
    # ONE Cholesky factorisation of the joint correlation matrix, with the
    # conditioning set ordered FIRST:
    #
    #   R = [[R11, R12], [R21, R22]] = L t(L),  L = [[L11, 0], [L21, L22]]
    #
    # Then with omega1 = L11^{-1} g1 and z ~ N(0, I),
    #
    #   g2 = L21 omega1 + L22 z
    #
    # is an exact draw from p(g2 | g1): its mean is L21 L11^{-1} g1 =
    # R21 R11^{-1} g1 and its covariance is L22 t(L22) =
    # R22 - R21 R11^{-1} R12. Neither R11^{-1} nor the conditional covariance
    # is ever formed.
    #
    # Everything here is on the correlation (omega) scale 
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

    # calculate discard probabilities and select "discarded" locations
    probs <- pnorm(-beta * omega_x_tilde_all)
    idx_keep <- rbinom(k_star, 1, probs) == 1

    disc_coords <- coords_pp[idx_keep, ]
    omega_x_tilde <- omega_x_tilde_all[idx_keep]

    all_coords <- rbind(obs_coords, disc_coords)
    return(list("all_coords" = all_coords, "omega_x_tilde" = omega_x_tilde))
  }

### Spatial Random Effects ####

#### Cholesky helpers ####

chol_field_factor <- function(coords, rho, jitter = 1e-6) {
  R <- exp_cor(fields::rdist(coords), range = rho) +
    diag(jitter, nrow(as.matrix(coords)))
  L <- FastGP::rcppeigen_get_chol(R)
  if (!all(is.finite(L))) {
    # fall back to base R if the Eigen factorisation fails
    L <- t(chol(R))
  }
  return(L)
}

# z = L^{-1} x, i.e. solves L z = x
whiten_chol <- function(L, x) {
  return(forwardsolve(L, as.numeric(x)))
}

# x' R^{-1} y (y defaults to x), from triangular solves only
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

### Spatial field + eta: joint elliptical slice sampler ####
#
# eta (one value per year) and omega(.) 
#
#   (eta - prior_eta_mean, omega(K)) ~ N(0, blockdiag(prior_eta_var * I_J, R(rho_S)))
#
# and q(s) = eta + sqrt(sigma2_S) * omega(s)
sample_eta_S_ess <- function(
    eta_cur,
    omega_n_cur,
    omega_x_tilde,
    beta,
    sigma2_S,
    rho_S,
    nu,
    xi,
    y,
    littlen,
    all_coords_new,
    stations,
    idx_years,
    prior_eta_mean = 0,
    prior_eta_var = 10,
    n_ess_passes = 20,
    chol_R = NULL) {
  k_new <- dim(all_coords_new)[1]
  n_years <- length(eta_cur)

  if (is.null(chol_R)) {
    chol_R <- chol_field_factor(all_coords_new, rho_S)
  }

  omega_k_cur <- c(omega_n_cur, omega_x_tilde)
  if (length(omega_k_cur) != k_new) {
    stop(
      "sample_eta_S_ess: field length does not match the number of locations"
    )
  }

  sd_eta <- sqrt(prior_eta_var)
  state_cur <- c(eta_cur - prior_eta_mean, omega_k_cur)

  lik_of <- function(state) {
    eta_vec <- state[1:n_years] + prior_eta_mean
    omega_vec <- state[(n_years + 1):(n_years + k_new)]
    lik_eta_omega(
      eta_vec,
      omega_vec,
      beta,
      sigma2_S,
      nu,
      xi,
      y,
      littlen,
      k_new,
      stations,
      idx_years
    )
  }

  lik_cur <- lik_of(state_cur)

  # Perform multiple passes of the ESS
  for (pass in seq_len(n_ess_passes)) {
    eps_eta <- rnorm(n_years, 0, sd_eta)
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

  eta_new <- state_cur[1:n_years] + prior_eta_mean
  omega_k_new <- state_cur[(n_years + 1):(n_years + k_new)]
  return(list(eta = eta_new, omega_k = omega_k_new))
}

# Joint log-target for (eta, omega): point-process term is sigma2_S-free
# under decomposing S (beta * omega, not (beta / sqrt(sigma2_S)) * S);
# the bGEV term uses q = eta + sqrt(sigma2_S) * omega but otherwise it's all the same
lik_eta_omega <-
  function(
      eta_vec,
      omega_k,
      beta,
      sigma2_S,
      nu,
      xi,
      y,
      littlen,
      k,
      stations,
      idx_years) {
    I_k <- c(rep(1, littlen), rep(-1, k - littlen))
    t1 <- sum(pnorm(beta * I_k * omega_k, log = TRUE))

    sigma_S <- sqrt(sigma2_S)
    q <- eta_vec[idx_years] + sigma_S * omega_k[stations]
    sb <- exp(nu[idx_years])
    xi_vec <- rep(xi, length(y))

    t2 <- sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
    return(t1 + t2)
  }

# bGEV part of the likelihood only (no point-process term). Used by the
# sigma2_S update, which never touches the point process under partial
# whitening.
lik_y_bgev <- function(S_vals, eta, nu, xi, y, stations, idx_years) {
  S_vals <- as.numeric(S_vals)
  q <- eta[idx_years] + S_vals[stations]
  sb <- exp(nu[idx_years])
  sum(dbgev2_vec(y, q, sb, rep(xi, length(y)), log = TRUE))
}

#### GEV Part ####

sample_nu <-
  function(
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
      iter,
      prev_prop_nu_var,
      default_var = 0.5,
      full_acc_rate,
      max_high = 3
    )
    nu_prop <- rnorm(1, nu_cur, prop_nu_var)
    log_mh_ratio <-
      target_nu(nu_prop, eta, xi, S_n, y, prior_nu_mean, prior_nu_var) -
      target_nu(nu_cur, eta, xi, S_n, y, prior_nu_mean, prior_nu_var)

    if (runif(1) < exp(log_mh_ratio)) {
      nu_acc <- 1
      nu_cur <- nu_prop
    }
    return(list(nu = nu_cur, nu_acc = nu_acc, prop_nu_var = prop_nu_var))
  }

target_nu <-
  function(nu_tmp, eta, xi, S_n, y, prior_nu_mean, prior_nu_var) {
    q <- eta + S_n
    sb <- rep(exp(nu_tmp), length(y))
    xi_vec <- rep(xi, length(y))

    t1 <- sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
    t2 <- (-1 / 2) * (nu_tmp - prior_nu_mean)^2 / prior_nu_var
    return(t1 + t2)
  }


sample_xi <-
  function(
      xi_cur,
      eta,
      nu,
      S_n,
      y,
      stations,
      idx_years,
      prior_xi_mean = -1,
      prior_xi_var = 1,
      iter = NULL,
      prev_prop_xi_var = NULL,
      full_acc_rate = NULL) {
    xi_acc <- 0
    prop_xi_var <- 0.7
    xi_prop <- rlnorm(1, log(xi_cur), prop_xi_var)
    num <-
      target_xi(
        xi_prop,
        eta,
        nu,
        S_n,
        y,
        stations,
        idx_years,
        prior_xi_mean,
        prior_xi_var
      )
    den <-
      target_xi(
        xi_cur,
        eta,
        nu,
        S_n,
        y,
        stations,
        idx_years,
        prior_xi_mean,
        prior_xi_var
      )
    log_mh_ratio <-
      num - den - log(xi_cur) + log(xi_prop)
    if (is.na(exp(log_mh_ratio))) {
      return(list(xi = xi_cur, xi_acc = xi_acc, prop_xi_var = prop_xi_var))
    }
    if (runif(1) < exp(log_mh_ratio)) {
      xi_acc <- 1
      xi_cur <- xi_prop
    }
    return(list(xi = xi_cur, xi_acc = xi_acc, prop_xi_var = prop_xi_var))
  }

target_xi <-
  function(
      xi_tmp,
      eta,
      nu,
      S_n,
      y,
      stations,
      idx_years,
      prior_xi_mean,
      prior_xi_var) {
    # Remember xi must always be positive.
    q <- eta[idx_years] + S_n[stations]
    sb <- exp(nu[idx_years])
    xi_vec <- rep(xi_tmp, length(y))

    t1 <- sum(dbgev2_vec(y, q, sb, xi_vec, log = TRUE))
    t2 <-
      dlnorm(xi_tmp, prior_xi_mean, sdlog = sqrt(prior_xi_var), log = TRUE)
    return(t1 + t2)
  }


#### Spatial Hyper-parameters Part ####

# sigma2_S : omega_k is held fixed, so neither the
# GP-density term nor the point-process term (both parameter-free in
# sigma2_S once omega is factored out) enter this ratio -- only the bGEV
# likelihood does.
sample_sigma2_S <-
  function(sigma2_S_cur,
           omega_k,
           eta,
           nu,
           xi,
           y,
           stations,
           idx_years,
           sigma2_S_prior_a,
           sigma2_S_prior_b,
           iter = NULL,
           prev_prop_sigma2_S_var = prev_prop_sigma2_S_var,
           full_acc_rate = sigma2_S_acc_rate) {
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
      lik_y_bgev(S_prop, eta, nu, xi, y, stations, idx_years) -
      lik_y_bgev(S_cur, eta, nu, xi, y, stations, idx_years) -
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
    # Proposal scale on the
    # log(rho) scale, adapted as for eta / nu / xi / beta.
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

#### Preferential Parameter Part (beta)####

# beta : only the point-process term ever involves
# beta, and that term is sigma2_S-free (beta * omega, not
# (beta / sqrt(sigma2_S)) * S), so sigma2_S no longer enters this step at all.
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
