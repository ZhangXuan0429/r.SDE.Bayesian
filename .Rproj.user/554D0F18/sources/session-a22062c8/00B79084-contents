#' @title Kernel ridge regression (KRR) with a user-supplied kernel
#'
#' @description Solve the linear system \eqn{(K + \lambda_{\mathrm{eff}} I)\alpha = y} and make predictions.
#'
#' @param X_train Numeric matrix (k x p). Training inputs.
#' @param y_train Numeric vector of length k. Training targets.
#' @param X_test Numeric matrix (m x p). Test inputs to predict on.
#' @param kernel_func Function taking (X1, X2, sigma) and returning a kernel matrix.
#' @param sigma Positive scalar bandwidth passed to `kernel_func`.
#' @param lambda Positive scalar base regularization.
#' @param gamma Positive scalar scaling for the ODE penalty; forms \eqn{\lambda_{\mathrm{eff}}}.
#'
#' @return A list with elements:
#' \describe{
#'   \item{alpha}{Numeric vector of coefficients \eqn{\alpha} (length k).}
#'   \item{y_pred}{Numeric vector of predictions on `X_test` (length m).}
#' }
#' @export
#'
#' @examples
#' X <- matrix(rnorm(20), nrow = 10); y <- rnorm(10)
#' fit <- krr(X, y, X, kernel_func = rbf_kernel, sigma = 1, lambda = 1e-2, gamma = 1)
#' head(fit$y_pred)
#'
krr <- function(X_train, y_train, X_test, kernel_func, sigma, lambda, gamma) {
  K_train <- kernel_func(X_train, X_train, sigma)
  n <- nrow(K_train)
  lambda_effective <- (lambda * n) / gamma

  K_reg <- K_train + lambda_effective * diag(n)

  # (K + λI)α = y_train
  alpha <- solve(K_reg, y_train, tol = 1e-10)

  K_test <- kernel_func(X_test, X_train, sigma)

  y_pred <- K_test %*% alpha

  return(list(alpha = alpha, y_pred = y_pred))
}
