#' Maximum Likelihood Estimation for Weibull Birnbaum-Saunders Distribution
#'
#' Estimates the parameters of the WBS distribution. This function is optimized
#' for simulation studies (fast execution, no console output, returns NAs on failure).
#'
#' @param y A numeric vector of observed data.
#' @return A numeric vector with the estimated parameters (alpha, beta, a, b).
#'         Returns a vector of NAs if optimization fails.
#' @examples
#'
#' alpha <- 1.5; beta <- 0.1; a <- 0.75; b <- 0.25
#' set.seed(123)
#' y <- rwbs(100, alpha, beta, a, b)
#' MLEwbs(y)
#' @export
MLEwbs<-function(y) {

  # 1. Optimized Initial Guess
  # Robust logic: beta close to the median, alpha close to standard deviation
  chute=log(c(alpha=mean(y),beta=sd(y),a=1,b=1) )

  # 2. Log-Likelihood Function
  # Note: Assumes dwbs() is available in the package namespace
  LKH <- function(par) {
    alpha <- exp(par[1])
    beta  <- exp(par[2])
    a     <- exp(par[3])
    b     <- exp(par[4])

    dens <- dwbs(y, alpha, beta, a, b)

    # Fast protection against invalid values (NaNs or Infinites) during optimization
    if (any(dens <= 0 | !is.finite(dens))) return(1e10)

    -sum(log(dens))
  }

  # 3. Optimization
  # Attempt 1: nlminb (Generally faster and more robust)
  resultado <- try(stats::nlminb(chute, LKH,
                                 lower = rep(-Inf, 4),
                                 upper = rep(Inf, 4)),
                   silent = TRUE)

  parametros <- NULL

  # Check for nlminb success
  if (!inherits(resultado, "try-error") && resultado$convergence == 0) {
    parametros <- resultado$par
    # Hessian calculation removed for speed in simulations
  } else {
    # Attempt 2: optim (BFGS) as fallback
    resultado <- try(stats::optim(
      par = chute,
      fn = LKH,
      method = "BFGS",
      control = list(maxit = 5000, reltol = 1e-12),
      hessian = FALSE # Hessian not needed for point estimation simulation
    ), silent = TRUE)

    if (inherits(resultado, "try-error")) {
      # In simulation, returning NA is better than stopping the script with an error
      return(c(alpha = NA, beta = NA, a = NA, b = NA))
    }

    parametros <- resultado$par
  }

  # 4. Transformation and Return
  MLE <- exp(parametros)
  names(MLE) <- c("alpha", "beta", "a", "b")

  return(MLE)
}
