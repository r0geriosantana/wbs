#' Weibull Birnbaum-Saunders Density Function
#'
#' @param t Time or failure time (t > 0)
#' @param alpha Shape parameter (alpha > 0)
#' @param beta Scale parameter (beta > 0)
#' @param a Shape parameter (a > 0)
#' @param b Shape parameter (b > 0)
#' @return Density value of Weibull-Birnbaum-Saunders distribution at t
#' @examples
#' dwbs(2.5,alpha=1.5,beta=3,a=2,b=5)
#' @export
#'
dwbs=function(t, alpha, beta, a, b) {
  # Função auxiliar para calcular v
  v <- function(t, alpha, beta) {
    alpha^(-1) * (sqrt(t / beta) - sqrt(beta / t))
  }

  # Constante k(alpha, beta) com exp(alpha^(-2))
  k <- function(alpha, beta) {
    exp(alpha^(-2)) / (2 * alpha * sqrt(2 * pi * beta))
  }

  # Function tau(z)
  tau <- function(z) {
    z + 1 / z
  }
  # Calcula v(t, alpha, beta)
  v.val <- v(t, alpha, beta)
  pv <- pnorm(v.val)  # Probabilidade normal acumulada

  # Parte 1: constante e termos iniciais
  part1=a * b * k(alpha, beta) * t^(-3/2) * (t + beta)
  # Parte 2: termos envolvendo pv
  part2=(pv^(b - 1)) / ((1 - pv)^(b + 1))
  # Parte 3: expoente final
  exponent= -tau(t / beta) / (2 * alpha^2) - a * (pv / (1 - pv))^b

  # Retorna a densidade
  part1 * part2 * exp(exponent)
}
