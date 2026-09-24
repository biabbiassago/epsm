# ===========================================================================
# Blended generalised extreme value (bGEV) distribution
# ===========================================================================
#
# Adapted from the R code accompanying:
#
#   Castro-Camilo, D., Huser, R. and Rue, H. (2022). Practical strategies for
#   generalized extreme value-based regression models for extremes.
#   Environmetrics, 33(6), e2742. https://doi.org/10.1002/env.2742
#   Code: https://github.com/dcastrocamilo/bGEV
#
# The bGEV distribution and its (q, s_beta, xi) parametrisation are from:
#
#   Vandeskog, S. M., Martino, S., Castro-Camilo, D. and Rue, H. (2022).
#   Modelling sub-daily precipitation extremes with the blended generalised
#   extreme value distribution. Journal of Agricultural, Biological and
#   Environmental Statistics, 27, 598-621.
#   https://doi.org/10.1007/s13253-022-00500-7
#
# Original functions by S. M. Vandeskog and D. Castro-Camilo. The *_vec
# functions (vectorised over observations, used by the MCMC samplers) are
# additions for this package.
#
# Parametrisations
#   "old": (mu, sigma, xi)  -- usual GEV location, scale, shape.
#   "new": (q, sb, xi)      -- q = alpha-quantile (median for alpha = 0.5),
#                              sb = spread q_{1-beta/2} - q_{beta/2}.
# The bGEV equals a Gumbel below the p_a quantile, the GEV above the p_b
# quantile, and a Beta(s, s)-weighted geometric mixture in between. Requires
# xi >= 0.
# ===========================================================================


# ---------------------------------------------------------------------------
# GEV building blocks (vectorised over all arguments)
# ---------------------------------------------------------------------------

# NOTE: all arguments are recycled to a common length first. ifelse()
# returns a result the length of its *test* (xi == 0), so in the original
# code a scalar xi (e.g. the 0 passed for the Gumbel component) with a
# vector x silently returned a single value, recycled.
pgev <- function(x, mu, sigma, xi) {
  n <- max(length(x), length(mu), length(sigma), length(xi))
  xi <- rep_len(xi, n)
  ifelse(
    xi == 0,
    exp(-exp(-(x - mu) / sigma)),
    exp(-pmax(0, 1 + xi * (x - mu) / sigma)^(-1 / xi))
  )
}

qgev <- function(p, mu, sigma, xi) {
  n <- max(length(p), length(mu), length(sigma), length(xi))
  xi <- rep_len(xi, n)
  ifelse(
    xi == 0,
    mu - sigma * log(-log(p)),
    mu - sigma * (1 / xi) * (1 - (-log(p))^(-xi))
  )
}

rgev <- function(n, mu, sigma, xi) {
  lengths <- sapply(list(mu, sigma, xi), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qgev(stats::runif(n), mu, sigma, xi)
}

# GEV density. Vectorised in xi (the original allbgev.R version used
# `if (xi == 0)` and so only accepted scalar xi). Accepts `log` (original
# name) or `log_scale` (name used by the *_vec functions).
dgev <- function(x, mu, sigma, xi, log = FALSE, log_scale = log) {
  n <- max(length(x), length(mu), length(sigma), length(xi))
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  xi <- rep_len(xi, n)

  z <- (x - mu) / sigma
  t <- pmax(1 + xi * z, 0)
  res <- ifelse(
    xi == 0,
    -exp(-z) - z,
    ifelse(t > 0, -t^(-1 / xi) - (1 / xi + 1) * log(t), -Inf)
  ) - log(sigma)
  if (!log_scale) res <- exp(res)
  res
}

# Return levels, GEV
return_level_gev <- function(period, mu, sigma, xi) {
  if (any(period <= 1)) warning("invalid period")
  p <- ifelse(period > 1, 1 - 1 / period, NA)
  qgev(p, mu, sigma, xi)
}

# Simulate N iid GEV (NOT bGEV) draws from (q, s, xi) parameters. Used by the
# simulation functions in sim-data.R. Formerly in reparametrized_gev.R, where
# rep2gev() did the conversion; it is identical to new.to.old().
# Uses SpatialExtremes::rgev when installed (so seeds reproduce the original  simulations), otherwise this file's rgev().
rgevrep <- function(N, q, s, xi, a = 0.5, b = 0.5) {
  p <- new.to.old(c(q, s, xi), alpha = a, beta = b)
  if (requireNamespace("SpatialExtremes", quietly = TRUE)) {
    SpatialExtremes::rgev(N, loc = p$mu, scale = p$sigma, shape = p$xi)
  } else {
    rgev(N, p$mu, p$sigma, p$xi)
  }
}

# GEV quantile function in the (q, s, xi) parametrisation.
qgevrep <- function(quant, q, s, xi, a = 0.5, b = 0.5) {
  p <- new.to.old(c(q, s, xi), alpha = a, beta = b)
  qgev(quant, p$mu, p$sigma, p$xi)
}

# GEV density in the (q, sb, xi) parametrisation (scalar parameters).
dgev2 <- function(x, q, sb, xi, alpha = 0.5, beta = 0.5, log = FALSE) {
  tmp <- new.to.old(c(q, sb, xi), alpha = alpha, beta = beta)
  dgev(x, tmp$mu, tmp$sigma, xi, log_scale = log)
}


# ---------------------------------------------------------------------------
# Parametrisation conversions
# ---------------------------------------------------------------------------

# (q, sb, xi) -> (mu, sigma, xi); `par` is a length-3 vector.
new.to.old <- function(par, alpha = 0.5, beta = 0.5) {
  q <- par[1]
  s <- par[2]
  xi <- par[3]
  if (xi == 0) {
    ell1 <- log(-log(alpha))
    ell2 <- log(-log(beta / 2))
    ell3 <- log(-log(1 - beta / 2))
    sigma <- s / (ell2 - ell3)
    mu <- q + sigma * ell1
  } else {
    ell1 <- (-log(alpha))^(-xi)
    ell2 <- (-log(beta / 2))^(-xi)
    ell3 <- (-log(1 - beta / 2))^(-xi)
    mu <- q - s * (ell1 - 1) / (ell3 - ell2)
    sigma <- xi * s / (ell3 - ell2)
  }
  list(mu = mu, sigma = sigma, xi = xi)
}

# Vectorised new.to.old: `par` is an n x 3 matrix with columns (q, sb, xi).
# Returns a list of length-n vectors mu, sigma, xi.
new_to_old_vec <- function(par, alpha = 0.5, beta = 0.5) {
  q <- par[, 1]
  s <- par[, 2]
  xi <- par[, 3]
  gum <- xi == 0

  ell1 <- ifelse(gum, log(-log(alpha)), (-log(alpha))^(-xi))
  ell2 <- ifelse(gum, log(-log(beta / 2)), (-log(beta / 2))^(-xi))
  ell3 <- ifelse(gum, log(-log(1 - beta / 2)), (-log(1 - beta / 2))^(-xi))

  sigma <- ifelse(gum, s / (ell2 - ell3), xi * s / (ell3 - ell2))
  mu <- ifelse(gum, q + sigma * ell1, q - s * (ell1 - 1) / (ell3 - ell2))
  list(mu = unname(mu), sigma = unname(sigma), xi = unname(xi))
}

# (mu, sigma, xi) -> (q, sb, xi); `par` is a length-3 vector.
old.to.new <- function(par, alpha = 0.5, beta = 0.5) {
  mu <- par[1]
  sigma <- par[2]
  xi <- par[3]
  qalpha <- qgev(alpha, mu, sigma, xi)
  qbeta1 <- qgev(beta / 2, mu, sigma, xi)
  qbeta2 <- qgev(1 - beta / 2, mu, sigma, xi)
  list(q = qalpha, s = qbeta2 - qbeta1, xi = xi)
}


# ---------------------------------------------------------------------------
# bGEV: distribution, quantile, density, random generation
# (original versions; arguments are recycled with fix_lengths())
# ---------------------------------------------------------------------------

pbgev <- function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  fix_lengths(x, mu, sigma, xi, p_a, p_b, s)
  g <- get_gumbel_par(mu, sigma, xi, p_a, p_b)
  a <- qgev(p_a, mu, sigma, xi)
  b <- qgev(p_b, mu, sigma, xi)
  p <- stats::pbeta((x - a) / (b - a), s, s)
  pgev(x, mu, sigma, xi)^p * pgev(x, g$mu, g$sigma, 0)^(1 - p)
}

qbgev <- function(p, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  fix_lengths(p, mu, sigma, xi, p_a, p_b, s)
  res <- rep(NA, length(p))
  gumbel <- which(p <= p_a)
  frechet <- which(p >= p_b)
  mixing <- which(p_a < p & p < p_b)
  if (any(gumbel)) {
    g <- get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] <- qgev(p[gumbel], g$mu, g$sigma, 0)
  }
  if (any(frechet)) {
    res[frechet] <- qgev(p[frechet], mu[frechet], sigma[frechet], xi[frechet])
  }
  if (any(mixing)) {
    res[mixing] <- qbgev_mixing(
      p[mixing], mu[mixing], sigma[mixing], xi[mixing],
      p_a[mixing], p_b[mixing], s[mixing]
    )
  }
  res
}

rbgev <- function(n, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  lengths <- sapply(list(mu, sigma, xi), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qbgev(stats::runif(n), mu, sigma, xi, p_a, p_b, s)
}

# Density, (mu, sigma, xi) parametrisation
dbgev <- function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  fix_lengths(x, mu, sigma, xi, p_a, p_b, s)
  a <- qgev(p_a, mu, sigma, xi)
  b <- qgev(p_b, mu, sigma, xi)
  res <- rep(NA, length(x))
  gumbel <- which(x <= a)
  frechet <- which(x >= b)
  mixing <- which(a < x & x < b)
  if (any(gumbel)) {
    g <- get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] <- dgev(x[gumbel], g$mu, g$sigma, 0, log_scale = log)
  }
  if (any(frechet)) {
    res[frechet] <- dgev(x[frechet], mu[frechet], sigma[frechet], xi[frechet], log_scale = log)
  }
  if (any(mixing)) {
    res[mixing] <- dbgev_mixing(
      x[mixing], mu[mixing], sigma[mixing], xi[mixing],
      p_a[mixing], p_b[mixing], s[mixing],
      log = log
    )
  }
  res
}

# Density, (q, sb, xi) parametrisation (scalar parameters)
dbgev2 <- function(x, q, sb, xi, alpha = 0.5, beta = 0.5,
                   p_a = .1, p_b = .2, s = 5, log = FALSE) {
  fix_lengths(x, q, sb, xi, p_a, p_b, s)
  tmp <- new.to.old(c(q, sb, xi), alpha = alpha, beta = beta)
  mu <- tmp$mu
  sigma <- tmp$sigma
  a <- qgev(p_a, mu, sigma, xi)
  b <- qgev(p_b, mu, sigma, xi)
  res <- rep(NA, length(x))
  gumbel <- which(x <= a)
  frechet <- which(x >= b)
  mixing <- which(a < x & x < b)
  if (any(gumbel)) {
    g <- get_gumbel_par(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] <- dgev(x[gumbel], g$mu, g$sigma, 0, log_scale = log)
  }
  if (any(frechet)) {
    res[frechet] <- dgev(x[frechet], mu[frechet], sigma[frechet], xi[frechet], log_scale = log)
  }
  if (any(mixing)) {
    res[mixing] <- dbgev_mixing(
      x[mixing], mu[mixing], sigma[mixing], xi[mixing],
      p_a[mixing], p_b[mixing], s[mixing],
      log = log
    )
  }
  res
}

# Return levels, (mu, sigma, xi) parametrisation
return_level_bgev <- function(period, mu, sigma, xi, p_a = .1, p_b = .2, s = 5) {
  if (any(period <= 1)) warning("invalid period")
  p <- ifelse(period > 1, 1 - 1 / period, NA)
  qbgev(p, mu, sigma, xi, p_a, p_b, s)
}

# Return levels, (q, sb, xi) parametrisation
return_level_bgev2 <- function(period, q, sb, xi, p_a = .1, p_b = .2, s = 5) {
  if (any(period <= 1)) warning("invalid period")
  p <- ifelse(period > 1, 1 - 1 / period, NA)
  old_pars <- new_to_old_vec(cbind(q, sb, xi))
  qbgev(p, old_pars$mu, old_pars$sigma, old_pars$xi, p_a, p_b, s)
}


# ---------------------------------------------------------------------------
# Vectorised bGEV (one parameter set per observation). dbgev2_vec is the
# likelihood used throughout the MCMC samplers.
# ---------------------------------------------------------------------------

dbgev2_vec <- function(x, q, sb, xi, alpha = 0.5, beta = 0.5,
                       p_a = rep(0.1, length(q)), p_b = rep(0.2, length(q)),
                       s = rep(5, length(q)), log = FALSE) {
  if (!(length(x) == length(q) && length(q) == length(sb) && length(sb) == length(xi))) {
    stop("x, q, sb and xi must all be vectors of the same length.")
  }

  tmp <- new_to_old_vec(cbind(q, sb, xi), alpha = alpha, beta = beta)
  mu <- tmp$mu
  sigma <- tmp$sigma

  a <- qgev(p_a, mu, sigma, xi)
  b <- qgev(p_b, mu, sigma, xi)

  res <- numeric(length(x))
  gumbel <- x <= a
  frechet <- x >= b
  mixing <- !gumbel & !frechet

  if (any(gumbel)) {
    g <- get_gumbel_par_vec(
      mu[gumbel], sigma[gumbel], xi[gumbel],
      p_a[gumbel], p_b[gumbel]
    )
    res[gumbel] <- dgev(x[gumbel], g$mu, g$sigma, rep(0, length(g$mu)), log_scale = log)
  }
  if (any(frechet)) {
    res[frechet] <- dgev(x[frechet], mu[frechet], sigma[frechet], xi[frechet], log_scale = log)
  }
  if (any(mixing)) {
    res[mixing] <- dbgev_mixing_vec(
      x[mixing], mu[mixing], sigma[mixing], xi[mixing],
      p_a[mixing], p_b[mixing], s[mixing],
      log_scale = log
    )
  }
  res
}

# Random generation, (q, sb, xi) parametrisation. n must be 1 (one uniform
# shared by all parameter sets) or length(q) (one draw per parameter set).
rbgev2 <- function(n, q, sb, xi) {
  if (!n %in% c(1, length(q))) stop("n must be 1 or length(q)")
  tmp <- new_to_old_vec(cbind(q, sb, xi))
  u <- if (n == 1) rep(stats::runif(1), length(q)) else stats::runif(n)
  qbgev_vec(u, tmp$mu, tmp$sigma, tmp$xi)
}

# Vectorised quantile function. Fixed relative to the original, which
# overwrote the whole result in each branch (so it was only correct when
# every p fell in the same region).
qbgev_vec <- function(p, mu, sigma, xi,
                      p_a = rep(.1, length(p)), p_b = rep(.2, length(p)),
                      s = rep(5, length(p))) {
  if (length(p) != length(mu)) stop("p and mu must be the same length")
  res <- rep(NA_real_, length(p))
  gumbel <- p <= p_a
  frechet <- p >= p_b
  mixing <- !gumbel & !frechet
  if (any(gumbel)) {
    g <- get_gumbel_par_vec(mu[gumbel], sigma[gumbel], xi[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] <- qgev(p[gumbel], g$mu, g$sigma, 0)
  }
  if (any(frechet)) {
    res[frechet] <- qgev(p[frechet], mu[frechet], sigma[frechet], xi[frechet])
  }
  if (any(mixing)) {
    res[mixing] <- qbgev_mixing_vec(
      p[mixing], mu[mixing], sigma[mixing], xi[mixing],
      p_a[mixing], p_b[mixing], s[mixing]
    )
  }
  res
}


# ---------------------------------------------------------------------------
# Helpers: Gumbel parameters, mixing-region density and quantile
# ---------------------------------------------------------------------------

# Parameters of the Gumbel G matching the GEV at the p_a and p_b quantiles.
get_gumbel_par <- function(mu, sigma, xi, p_a = .1, p_b = .2) {
  if (any(xi < 0)) stop("xi must be nonnegative")
  a <- qgev(p_a, mu, sigma, xi)
  b <- qgev(p_b, mu, sigma, xi)
  sigma2 <- (b - a) / log(log(p_a) / log(p_b))
  mu2 <- a + sigma2 * log(-log(p_a))
  list(mu = mu2, sigma = sigma2)
}

get_gumbel_par_vec <- function(mu, sigma, xi,
                               p_a = rep(.1, length(mu)), p_b = rep(.2, length(mu))) {
  get_gumbel_par(mu, sigma, xi, p_a, p_b)
}

# Density on the mixing region a < x < b
dbgev_mixing <- function(x, mu, sigma, xi, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  g <- get_gumbel_par(mu, sigma, xi, p_a, p_b)
  a <- qgev(p_a, mu, sigma, xi)
  b <- qgev(p_b, mu, sigma, xi)
  if (any(x <= a | x >= b)) stop("x is outside the domain for mixing")
  p <- stats::pbeta((x - a) / (b - a), s, s)
  p_der <- stats::dbeta((x - a) / (b - a), s, s) / (b - a)
  term1 <- -p_der * (1 + xi * (x - mu) / sigma)^(-1 / xi)
  term2 <- p / sigma * (1 + xi * (x - mu) / sigma)^(-1 / xi - 1)
  term3 <- p_der * exp(-(x - g$mu) / g$sigma)
  term4 <- (1 - p) / g$sigma * exp(-(x - g$mu) / g$sigma)
  term0 <- p * log(pgev(x, mu, sigma, xi)) + (1 - p) * log(pgev(x, g$mu, g$sigma, 0))
  res <- term0 + log(term1 + term2 + term3 + term4)
  if (!log) res <- exp(res)
  res
}

dbgev_mixing_vec <- function(x, mu, sigma, xi,
                             p_a = rep(.1, length(mu)), p_b = rep(.2, length(mu)),
                             s = rep(5, length(mu)), log_scale = FALSE) {
  dbgev_mixing(x, mu, sigma, xi, p_a, p_b, s, log = log_scale)
}

# Quantile on the mixing region p_a < p < p_b (numerical inversion)
qbgev_mixing <- function(p, mu, sigma, xi, p_a = .1, p_b = .2, s = 5,
                         lower = 0, upper = 100) {
  if (any(p <= p_a | p >= p_b)) stop("p is outside the domain for mixing")
  res <- vector("numeric", length(p))
  for (i in seq_along(p)) {
    f <- function(x) (pbgev(x, mu, sigma, xi, p_a, p_b, s) - p)[i]
    sol <- stats::uniroot(f, lower = lower, upper = upper, extendInt = "upX")
    res[i] <- sol$root
  }
  res
}

qbgev_mixing_vec <- function(p, mu, sigma, xi,
                             p_a = rep(.1, length(p)), p_b = rep(.2, length(p)),
                             s = rep(5, length(p)), lower = 0, upper = 100) {
  if (any(p <= p_a | p >= p_b)) stop("p is outside the domain for mixing")
  res <- vector("numeric", length(p))
  for (i in seq_along(p)) {
    f <- function(x) pbgev(x, mu[i], sigma[i], xi[i], p_a[i], p_b[i], s[i]) - p[i]
    sol <- stats::uniroot(f, lower = lower, upper = upper, extendInt = "upX")
    res[i] <- sol$root
  }
  res
}

# Recycle the named arguments to a common length, in the caller's frame.
fix_lengths <- function(...) {
  call <- match.call()
  varnames <- sapply(call[-1], as.character)
  e <- parent.frame()
  vars <- lapply(varnames, get, envir = e)
  lengths <- sapply(vars, length)
  max_length <- max(lengths)
  if (any(max_length %% lengths != 0)) stop("Bad input lengths")
  for (i in seq_along(vars)) {
    if (lengths[i] < max_length) {
      assign(varnames[i], rep(vars[[i]], max_length / lengths[i]), envir = e)
    }
  }
  0
}


# ---------------------------------------------------------------------------
# Backward-compatible names. These duplicate functions above but are kept so
# existing scripts that call them keep working.
# ---------------------------------------------------------------------------

# reparametrized_gev.R
elle <- function(a, xi) {
  if (xi != 0) (-log(a))^(-xi) else log(-log(a))
}

rep2gev <- function(q, s, xi, a = 0.5, b = 0.5) {
  new.to.old(c(q, s, xi), alpha = a, beta = b)
}

# (Original had an undefined `x` in the xi == 0 branch; fixed.)
gev2rep <- function(mu, sigma, xi, a = 0.5, b = 0.5) {
  if (xi != 0) {
    den <- elle(1 - b / 2, xi) - elle(b / 2, xi)
    q <- mu + sigma * (elle(a, xi) - 1) / xi
    s <- sigma * den / xi
  } else {
    den <- elle(b / 2, xi) - elle(1 - b / 2, xi)
    q <- mu - sigma * elle(a, xi)
    s <- sigma * den
  }
  list(q = q, s = s, xi = xi)
}

# GEVcode.R
dgev_base <- function(x, mu, sigma, xi, log = FALSE) {
  dgev(x, mu, sigma, xi, log_scale = log)
}

qgev_other <- function(p, mu, sigma, xi) qgev(p, mu, sigma, xi)

# bGEVcode.R
dbgev2_int <- function(x, q, sb, xi, alpha = 0.5, beta = 0.5,
                       p_a = 0.1, p_b = 0.2, s = 5, log = FALSE) {
  dbgev2(x, q, sb, xi, alpha, beta, p_a, p_b, s, log)
}
