# Check simplex constraint x1 + x2 = 1
.sde_check_simplex <- function(x, tol = 1e-8) {
  x <- as.numeric(x)
  if (length(x) != 2) return(FALSE)
  if (any(!is.finite(x))) return(FALSE)
  if (any(x < -tol) || any(x > 1 + tol)) return(FALSE)
  if (abs(sum(x) - 1) > tol) return(FALSE)
  TRUE
}

#' Log conditional density under Euler--Maruyama Gaussian approximation
#'
#' Computes \eqn{\log p(x_{t+\Delta t} \mid x_t, \theta)} under an Euler--Maruyama
#' Gaussian approximation for a two-sense LV-style competition SDE.
#'
#' The model is written as
#' \deqn{
#' dx_1 = b_1(x)\,dt + c\,dW,\quad
#' dx_2 = b_2(x)\,dt - c\,dW,
#' }
#' and we enforce the simplex constraint \eqn{x_1 + x_2 = 1}.
#'
#' Because the state is constrained on the simplex (x1 + x2 = 1), we compute the
#' Euler--Maruyama transition density for x1 and verify that x2 = 1 - x1 holds
#' (up to tolerance).
#'
#' The approximation is:
#' \deqn{
#' x_{1,t+dt} | x_{1,t} \sim N(x_{1,t} + b_1(x_t) dt,\; dt c^2).
#' }

#'
#' @param xt1 Numeric length-2 vector, state at time t: c(x1, x2).
#' @param xt2 Numeric length-2 vector, state at time t+dt: c(x1, x2).
#' @param dt  Positive numeric scalar, time step.
#' @param drift Numeric length-2 vector: c(\eqn{(r_1, \sigma_1)}).
#' @param diffusion Positive numeric scalar: c.
#' @param tol Nonnegative numeric scalar, tolerance for simplex checks.
#' @return Numeric scalar, log conditional density.
#' @export
sde_log_cond_pdf <- function(xt1, xt2, dt, drift, diffusion, tol = 1e-8) {
  xt1 <- as.numeric(xt1)
  xt2 <- as.numeric(xt2)
  if (length(xt1) != 2 || length(xt2) != 2) stop("xt1/xt2 must be length-2 vectors.")
  if (!is.finite(dt) || dt <= 0) return(-1e10)

  r1 <- as.numeric(drift[1])
  s1 <- as.numeric(drift[2])
  c  <- as.numeric(diffusion)

  if (!all(is.finite(c(r1, s1, c)))) return(-1e10)
  if (r1 < 0 || s1 < 0) return(-1e10)
  if (c <= 0) return(-1e10)

  # Enforce simplex for proportions
  if (!.sde_check_simplex(xt1, tol = tol)) return(-1e10)
  if (!.sde_check_simplex(xt2, tol = tol)) return(-1e10)
  if (abs(xt2[2] - (1 - xt2[1])) > 10 * tol) return(-1e10)

  x1 <- xt1[1]
  x2 <- xt1[2]

  # LV-style drift for x1: b1 = r1*(1-sigma1)*x1*x2
  b1 <- r1 * (1 - s1) * x1 * x2

  mu1  <- x1 + b1 * dt
  var1 <- dt * c * c
  if (!is.finite(var1) || var1 <= 0) return(-1e10)

  z <- xt2[1] - mu1
  logdens <- -0.5 * log(2 * pi * var1) - 0.5 * (z * z) / var1
  as.numeric(logdens)
}

#' Log posterior for SDE parameters (Uniform prior on drift only)
#'
#' Parameter vector is \eqn{(r_1, \sigma_1, c)}.
#' A Uniform(0, prior_max) prior is applied to \eqn{(r_1, \sigma_1)} only.
#' No upper-bound prior is imposed on \eqn{c} (we only enforce \eqn{c > 0}).
#'
#' The time step \eqn{\Delta t} is computed from the time column in \code{dat}.
#'
#' @param par Numeric vector length 3: c\eqn{(r_1, \sigma_1, c)}.
#' @param dat Data frame or matrix with columns: time, x1, x2.
#' @param prior_max Positive numeric scalar, upper bound for Uniform prior on drift.
#' @param log_cond_pdf Function computing log conditional density.
#' @param tol Tolerance for simplex checks.
#'
#' @return Numeric scalar, log posterior.
#' @export
sde_log_posterior <- function(par,
                              dat,
                              prior_max = 100,
                              log_cond_pdf = sde_log_cond_pdf,
                              tol = 1e-8) {
  par <- as.numeric(par)
  if (length(par) != 3) stop("par must be length 3: c(r1, sigma1, c).")
  if (!is.finite(prior_max) || prior_max <= 0) stop("prior_max must be positive.")

  drift <- as.numeric(par[1:2])
  c     <- as.numeric(par[3])

  # Uniform prior for drift only
  lp_drift <- sum(stats::dunif(drift, min = 0, max = prior_max, log = TRUE))
  if (!is.finite(lp_drift)) return(-1e10)
  if (!is.finite(c) || c <= 0) return(-1e10)

  dat <- as.data.frame(dat)
  if (ncol(dat) < 3) stop("dat must have columns: time, x1, x2.")

  tvec <- as.numeric(dat[[1]])
  X <- as.matrix(dat[, 2:3, drop = FALSE])

  if (nrow(X) < 2) stop("dat must contain at least two time points.")
  if (any(!is.finite(tvec)) || any(!is.finite(X))) return(-1e10)

  dt_vec <- diff(tvec)
  if (any(!is.finite(dt_vec)) || any(dt_vec <= 0)) return(-1e10)

  ll <- 0
  for (i in seq_len(nrow(X) - 1)) {
    ll <- ll + log_cond_pdf(
      xt1 = X[i, ],
      xt2 = X[i + 1, ],
      dt  = dt_vec[i],
      drift = drift,
      diffusion = c,
      tol = tol
    )
  }

  as.numeric(ll + lp_drift)
}


#' Fit SDE parameters via MCMC (fmcmc)
#'
#' Runs MCMC for \eqn{(r_1, \sigma_1, c)} using \code{sde_log_posterior()}.
#' You may supply a different \code{log_cond_pdf} to fit alternative drift forms
#' (e.g., replicator dynamics) under the same inference pipeline.
#'
#' @param dat Data frame/matrix with columns: time, x1, x2.
#' @param initial Numeric vector length 3, initial parameters.
#' @param nsteps Integer, number of MCMC iterations.
#' @param prior_max Upper bound for Uniform prior on drift parameters.
#' @param log_cond_pdf Function computing log conditional density.
#' @param tol Tolerance for simplex checks.
#' @param kernel Proposal kernel from fmcmc, e.g., fmcmc::kernel_ram().
#'
#' @return MCMC chain matrix from fmcmc::MCMC().
#' @export
#' @examples
#' # Example data: one polysemous word with two-sense proportion trajectory
#' result <- polysemous_data$entertain
#' dat <- result$df  # columns: time, x1, x2
#'
#' set.seed(123)
#' # MCMC can be time-consuming; reduce steps for a quick example.
#' chain <- sde_mcmc_fit(
#'   dat = dat,
#'   initial = rep(0.01, 3),
#'   nsteps = 2000,
#'   prior_max = 100,
#'   log_cond_pdf = sde_log_cond_pdf,
#'   tol = 1e-8,
#'   kernel = fmcmc::kernel_ram()
#' )
#'
#'# Log-posterior trace
#' lp <- apply(chain, 1, function(p) sde_log_posterior(p, dat = dat, prior_max = 100))
#'
#' # Posterior mean of the last samples
#' par_mean <- colMeans(tail(chain, 500))
#' par_mean
sde_mcmc_fit <- function(dat,
                         initial = rep(0.01, 3),
                         nsteps = 20000,
                         prior_max = 100,
                         log_cond_pdf = sde_log_cond_pdf,
                         tol = 1e-8,
                         kernel = fmcmc::kernel_ram()) {
  if (!requireNamespace("fmcmc", quietly = TRUE)) {
    stop("Package 'fmcmc' is required. Please install it.")
  }

  target <- function(par) {
    sde_log_posterior(
      par = par,
      dat = dat,
      prior_max = prior_max,
      log_cond_pdf = log_cond_pdf,
      tol = tol
    )
  }

  fmcmc::MCMC(fun = target, initial = initial, nsteps = as.integer(nsteps), kernel = kernel)
}
