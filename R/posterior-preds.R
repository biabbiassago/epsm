### """Predict S and y at unmeasured locations""" ###
# source(here::here("R/bGEV/utils.R"))

# call once, outside the MCMC loop (only if coords_pred is fixed)
precompute_pred_dists <- function(coords_pred) {
  fields::rdist(coords_pred)
}

# call every iteration
get_S_predlocs <- function(
    rho_S,
    sigma2_S,
    S_k_cur,
    coords_pred,
    coords_cur,
    dist_11 = NULL) {
  if (is.null(dist_11)) {
    dist_11 <- fields::rdist(coords_pred)
  }

  dist_12 <- fields::rdist(coords_pred, coords_cur)
  dist_22 <- fields::rdist(coords_cur)

  R_11 <- exp_cor(dist_11, rho_S)
  R_12 <- exp_cor(dist_12, rho_S)
  R_22 <- exp_cor(dist_22, rho_S)

  U <- chol(R_22)
  rhs <- cbind(S_k_cur, t(R_12))
  X <- backsolve(U, forwardsolve(t(U), rhs))

  cond_mean <- R_12 %*% X[, 1, drop = FALSE]
  cond_var <- sigma2_S * (R_11 - R_12 %*% X[, -1, drop = FALSE])

  FastGP::rcpp_rmvnorm(1, cond_var, cond_mean)
}

#
# get_fast_S_predlocs <- function(
#     rho_S,
#     sigma2_S,
#     S_k_cur,
#     coords_pred,
#     coords_cur,
#     dist_type) {
#   k_pred <- dim(coords_pred)[1]
#   k_cur <- dim(coords_cur)[1]
#
#   # Compute distance matrix in R
#   dist_big_coords <- if (dist_type == "euclidean") {
#     fields::rdist(rbind(coords_pred, coords_cur))
#   } else {
#     fields::rdist.earth(rbind(coords_pred, coords_cur))
#   }
#
#   # Compute correlation matrix
#   big_R <- exp_cor(dist_big_coords, rho_S)
#
#   # Partition the correlation matrix
#   R_11 <- big_R[1:k_pred, 1:k_pred]
#   R_12 <- big_R[1:k_pred, (k_pred + 1):(k_pred + k_cur)]
#   R_21 <- t(R_12)
#   R_22 <- big_R[(k_pred + 1):(k_pred + k_cur), (k_pred + 1):(k_pred + k_cur)]
#
#   # Call Rcpp function for expensive matrix computations & MVN sampling
#   S_pred_iter <- fast_mvn_conditional(R_11, R_12, R_21, R_22, sigma2_S, S_k_cur)
#
#   return(t(S_pred_iter))
# }

get_T_predlocs <- function(
    rho_T,
    sigma2_T,
    T_n,
    coords_pred,
    coords_obs,
    dist_type) {
  k_pred <- dim(coords_pred)[1]
  k_cur <- dim(coords_obs)[1]
  dist_big_coords <- if (dist_type == "euclidian") {
    fields::rdist(rbind(coords_pred, coords_obs))
  } else {
    fields::rdist.earth(rbind(coords_pred, coords_obs))
  }
  big_R <-
    exp_cor(dist_big_coords, rho_T)
  R_11 <- big_R[1:k_pred, 1:k_pred]
  R_12 <- big_R[1:k_pred, (k_pred + 1):(k_pred + k_cur)]
  R_21 <- t(R_12)
  R_22 <- big_R[(k_pred + 1):(k_pred + k_cur), (k_pred + 1):(k_pred + k_cur)]
  R_22_inv <- chol2inv(chol(R_22))
  cond_mean <- R_12 %*% R_22_inv %*% T_n
  cond_var <- sigma2_T * (R_11 - R_12 %*% R_22_inv %*% R_21)

  T_pred_iter <- mvrnorm(1, cond_mean, cond_var)

  return(T_pred_iter)
}
get_y_predlocs <- function(eta, S_pred, nu, T_pred = NULL, xi) {
  # to do after the run has finished
  k_pred <- dim(S_pred)[1]
  y_pred <- matrix(NA, nrow = k_pred, ncol = length(eta))
  if (is.null(T_pred)) {
    for (j in 1:k_pred) {
      pars <- matrix(
        c("q" = eta + S_pred[j, ], "s" = exp(nu), "xi" = xi),
        nrow = 3,
        ncol = length(eta),
        byrow = T
      )
      old.pars <- apply(
        pars,
        2,
        function(x) {
          new.to.old(c(x[1], x[2], x[3]))
        },
        simplify = T
      )
      y_pred[j, ] <- unlist(lapply(old.pars, function(x) {
        rbgev(1, mu = x$mu, sigma = x$sigma, xi = x$xi)
      }))
    }
  } else {
    for (j in 1:k_pred) {
      pars <- matrix(
        c("q" = eta + S_pred[j, ], "s" = exp(nu + T_pred[j, ]), "xi" = xi),
        nrow = 3,
        ncol = length(eta),
        byrow = T
      )
      old.pars <- apply(
        pars,
        2,
        function(x) {
          new.to.old(c(x[1], x[2], x[3]))
        },
        simplify = T
      )
      y_pred[j, ] <- unlist(lapply(old.pars, function(x) {
        rbgev(1, mu = x$mu, sigma = x$sigma, xi = x$xi)
      }))
    }
  }
  return(y_pred)
}
