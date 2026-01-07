#' Replicator SDE: log conditional density (Euler--Maruyama)
#'
#' Two-sense replicator SDE with shared parameters:
#' \deqn{
#' dx_1 = x_1 (f_1 - \bar f)\, dt + c\, dW,\quad
#' dx_2 = x_2 (f_2 - \bar f)\, dt - c\, dW,
#' }
#' where
#' \deqn{
#' f_1 = r (1 - x_1 - s x_2),\quad
#' f_2 = r (1 - x_2 - s x_1),\quad
#' \bar f = x_1 f_1 + x_2 f_2.
#' }
#'
#' Because x1 + x2 = 1, we compute the transition density for x1 and check the
#' simplex constraint:
#' \deqn{
#' x_{1,t+dt} | x_{1,t} \sim N(x_{1,t} + b_1(x_t) dt,\; dt c^2),
#' }
#' with \eqn{b_1(x_t) = x_1 (f_1 - \bar{f})}.
#'
#' @param xt1 Numeric length-2 vector: c(x1, x2).
#' @param xt2 Numeric length-2 vector: c(x1, x2).
#' @param dt Positive numeric scalar, time step.
#' @param drift Numeric length-2 vector: c(r, s).
#' @param diffusion Positive numeric scalar: c.
#' @param tol Nonnegative numeric scalar, tolerance for simplex checks.
#' @return Numeric scalar, log conditional density.
#' @export
#' @examples
#' # English comments only
#' result <- it1ws5_data_v2$entertain
#' dat <- result$df
#'
#' set.seed(123)
#' chain <- sde_mcmc_fit(
#'   dat = dat,
#'   initial = rep(0.01, 3),
#'   nsteps = 2000,
#'   prior_max = 100,
#'   log_cond_pdf = sde_log_cond_pdf_replica
#' )
#'
#' lp <- apply(chain, 1, function(p) sde_log_posterior(
#'   par = p, dat = dat, prior_max = 100, log_cond_pdf = sde_log_cond_pdf_replica
#' ))
#' par_mean <- colMeans(tail(chain, 500))
#' par_mean

sde_log_cond_pdf_replica <- function(xt1, xt2, dt, drift, diffusion, tol = 1e-8) {
  xt1 <- as.numeric(xt1)
  xt2 <- as.numeric(xt2)
  if (length(xt1) != 2 || length(xt2) != 2) stop("xt1/xt2 must be length-2 vectors.")
  if (!is.finite(dt) || dt <= 0) return(-1e10)

  r <- as.numeric(drift[1])
  s <- as.numeric(drift[2])
  c <- as.numeric(diffusion)

  if (!all(is.finite(c(r, s, c)))) return(-1e10)
  if (r < 0 || s < 0) return(-1e10)
  if (c <= 0) return(-1e10)

  if (!.sde_check_simplex(xt1, tol = tol)) return(-1e10)
  if (!.sde_check_simplex(xt2, tol = tol)) return(-1e10)
  if (abs(xt2[2] - (1 - xt2[1])) > 10 * tol) return(-1e10)

  x1 <- xt1[1]
  x2 <- xt1[2]

  f1 <- r * (1 - x1 - s * x2)
  f2 <- r * (1 - x2 - s * x1)
  fbar <- x1 * f1 + x2 * f2

  b1 <- x1 * (f1 - fbar)

  mu1  <- x1 + b1 * dt
  var1 <- dt * c * c
  if (!is.finite(var1) || var1 <= 0) return(-1e10)

  z <- xt2[1] - mu1
  logdens <- -0.5 * log(2 * pi * var1) - 0.5 * (z * z) / var1
  as.numeric(logdens)
}


