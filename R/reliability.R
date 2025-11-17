#' Calculate the WBS Reliability Function
#'
#' Computes the reliability function (survival function) R(t) for the
#' Weibull Birnbaum-Saunders (WBS) distribution.
#'
#' @param t A vector of positive quantiles (time).
#' @param alpha The shape parameter 'alpha' (must be positive).
#' @param beta The scale parameter 'beta' (must be positive).
#' @param a The shape parameter 'a' (must be positive).
#' @param b The shape parameter 'b' (must be positive).
#'
#' @return A vector of probabilities, R(t) = P(T > t).
#'
#' @examples
#' reliability(t = 2.5, alpha = 1.5, beta = 0.1, a = 0.75, b = 0.25)
#' @export
#'
reliability <- function(t, alpha, beta, a, b) {
  v = function(t, alpha, beta) {
    p = function(z) { z^(1/2) - z^(-1/2) }
    alpha^(-1) * p(t / beta)
  }
  num <- pnorm(v(t, alpha, beta), 0, 1)
  den <- 1 - pnorm(v(t, alpha, beta), 0, 1)
  # Reliability calculation
  exp(-a * ((num / den)^b))
}
