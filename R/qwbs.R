#' Weibull Birnbaum-Saunders Quantile Function
#'
#' @param u Probability value (u in [0, 1])
#' @param alpha Shape parameter (alpha > 0)
#' @param beta Scale parameter (beta > 0)
#' @param a Shape parameter (a > 0)
#' @param b Shape parameter (b > 0)
#' @return Quantile corresponding to probability u
#' @examples
#' qwbs(0.5, 0.5, 1, 1, 2)
#' @export
#'
qwbs=function(u, alpha, beta, a, b) {
  if (any(c(u) < 0 | u > 1)) {
    stop("u must be between 0 and 1.")
  }

  p <- (-log(1 - u) / a)^(1 / b) / (1 + (-log(1 - u) / a)^(1 / b))
  z <- qnorm(p)
  quantile <- (beta / 2) * (alpha * z + sqrt(4 + (alpha * z)^2))^2
  return(quantile)
}
