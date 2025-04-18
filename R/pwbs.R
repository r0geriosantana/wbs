#' Weibull Birnbaum-Saunders Cumulative Distribution Function
#'
#' @param t Time or failure time (t > 0)
#' @param alpha Shape parameter (alpha > 0)
#' @param beta Scale parameter (beta > 0)
#' @param a Shape parameter (a > 0)
#' @param b Shape parameter (b > 0)
#' @return CDF value of Weibull Birnbaum-Saunders distribution at t
#' @examples
#' pwbs(0.5,alpha=1.5,beta=2,a=0.5,b=1)
#' @export
#'
pwbs=function(t, alpha, beta, a, b) {
  if (any(c(t, alpha, beta, a, b) <= 0)) {
    stop("All parameters must be positive.")
  }

  # v function for transformation
  v <- function(t, alpha, beta) {
    alpha^(-1) * (sqrt(t / beta) - sqrt(beta / t))
  }

  v_val <- v(t, alpha, beta)
  q_v <- pnorm(v_val)  # CDF of normal at v

  # CDF formula for Weibull Birnbaum-Saunders
  1 - exp(-a * (q_v / (1 - q_v))^b)
}
