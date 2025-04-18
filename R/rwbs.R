#' Random Number Generation for Weibull Birnbaum-Saunders Distribution
#'
#' @param n Number of random values to generate
#' @param alpha Shape parameter (alpha > 0)
#' @param beta Scale parameter (beta > 0)
#' @param a Shape parameter (a > 0)
#' @param b Shape parameter (b > 0)
#' @return A vector of n random values following the WBS distribution
#' @examples
#' rwbs(10, 0.5, 1, 1, 2)
#'
#' @export
#'
rwbs=function(n, alpha, beta, a, b) {
  if (any(c(alpha, beta, a, b) <= 0)) {
    stop("All parameters must be positive.")
  }
  u=runif(n)
  quantiles=qwbs(u, alpha, beta, a, b)
  return(quantiles)
}

