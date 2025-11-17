#' Calculate the WBS Hazard Rate Function
#'
#' Computes the hazard rate function (failure rate) h(t) for the
#' Weibull Birnbaum-Saunders (WBS) distribution.
#'
#' @param t A vector of positive quantiles (time).
#' @param alpha The shape parameter 'alpha' (must be positive).
#' @param beta The scale parameter 'beta' (must be positive).
#' @param a The shape parameter 'a' (must be positive).
#' @param b The shape parameter 'b' (must be positive).
#'
#' @return A vector of hazard rates, h(t).
#'
#' @examples
#' failure(t = 2.5, alpha = 1.5, beta = 0.1, a = 0.75, b = 0.25)
#' @export
#'
failure <- function(t, alpha, beta, a, b) {
  # pt1: Constants and initial terms
  v = function(t, alpha, beta) {
    p = function(z) { z^(1/2) - z^(-1/2) }
    alpha^(-1) * p(t / beta)
  }
  # pt2: Kappa function
  k = function(alpha, beta) {
    exp(alpha^(-2)) / (2 * alpha * sqrt(2 * pi * beta))
  }
  # pt3: Tau function
  tau = function(z) { z + 1 / z }

  # pt4: Main term part A
  pt5 = a * b * k(alpha, beta) * (t^(-3/2)) * (t + beta)

  # pt5: Phi(v) terms
  v_val = v(t, alpha, beta)
  q = function(val) { pnorm(val, 0, 1) }

  # pt6: CDF ratio term
  pt6 = (q(v_val)^(b - 1)) / ((1 - q(v_val))^(b + 1))

  # pt7: Combine A and B
  pt7 = pt5 * pt6

  # pt8: Exponential term
  pt8 = -tau(t / beta) / (2 * alpha^2)

  # Final result
  pt7 * exp(pt8)
}
