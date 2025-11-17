#' Maximum Likelihood Estimation for Weibull Birnbaum-Saunders Distribution
#'
#' @param y A numeric vector of observed data.
#' @param plot Logical; if TRUE, plots the Q-Q plot of residuals.
#' @return A list containing the estimated parameters, AIC, BIC, and a Q-Q plot of residuals.
#' @examples
#'
#' alpha=1.5; beta=0.1; a=0.75;b=0.25
#' set.seed(123)
#' y=rwbs(100,alpha,beta,a,b)
#' modelwbs(y,plot='FALSE')
#' @export
#'
modelwbs=function(y,plot='FALSE'){
  chute=log(c(alpha=mean(y),beta=sd(y),a=1,b=1) )
  # Log-Likelihood (Negative) Function
  # Note: Assumes dwbs() is available in the package namespace
  LKH <- function(par) {
    alpha <- exp(par[1])
    beta  <- exp(par[2])
    a     <- exp(par[3])
    b     <- exp(par[4])

    # Protection against NaNs or Infinites during optimization.
    dens <- dwbs(y, alpha, beta, a, b)
    if (any(dens <= 0 | !is.finite(dens))) return(1e10)

    -sum(log(dens))
  }

  # 1. Optimization
  # Try nlminb first
  resultado <- try(stats::nlminb(chute, LKH,
                                 lower = rep(-Inf, 4),
                                 upper = rep(Inf, 4)),
                   silent = TRUE)

  # Check nlminb for success
  usar_optim <- inherits(resultado, "try-error") || resultado$convergence != 0

  if (!usar_optim) {
    message("Successful optimization with 'nlminb'.")
    parametros <- resultado$par
    log_vero   <- -resultado$objective # nlminb minimizes the negative, we invert it to LogLik

    # Calculate the Hessian function manually (since nlminb does not return it).
    hessiana <- numDeriv::hessian(func = LKH, x = parametros)

  } else {
    message("'nlminb' failed or did not converge. Trying with 'optim' (BFGS)...")

    resultado <- try(stats::optim(
      par = chute,
      fn = LKH,
      method = "BFGS",
      control = list(maxit = 5000, reltol = 1e-12),
      hessian = TRUE
    ), silent = TRUE)

    if (inherits(resultado, "try-error")) {
      stop("Optimization failed with both methods.")
    }

    parametros <- resultado$par
    log_vero   <- -resultado$value # optim minimizes the negative
    hessiana   <- resultado$hessian
  }

  # 2. Calculation of the Covariance Matrix and Standard Errors
  # Attempts to invert the Hessian
  Vdelta <- tryCatch({
    solve(hessiana)
  }, error = function(e) {
    warning("Non-invertible Hessian. Standard errors not calculated.")
    return(NULL)
  })

  # 3. AIC and BIC
  posto <- length(parametros)
  aic_val <- 2 * -log_vero + 2 * posto # Note teh sign: AIC = 2k - 2ln(L)
  bic_val <- 2 * -log_vero + posto * log(length(y))

  tabela_info <- round(data.frame(
    AIC = aic_val,
    BIC = bic_val,
    Log_likelihood = log_vero,
    N_Obs = length(y)
  ), 3)

  cat("\n--- Model Information ---\n")
  print(tabela_info)

  # 4. Confidence Interval (Delta Method)
  tabela_ic <- NULL
  estimativas_originais <- NULL

  # Transformation function g(even) = exp(even)
  g_transform <- function(par) {
    exp(par)
  }

  estimativas_originais <- g_transform(parametros)

  if (!is.null(Vdelta)) {
    # Jacobiana of transformation
    G_jacobiana <- numDeriv::jacobian(func = g_transform, x = parametros)

    # Variance in the original scale: J * V * J'
    V_p_mu <- G_jacobiana %*% Vdelta %*% t(G_jacobiana)

    # Standard errors on the original scale
    erros_padrao_orig <- sqrt(diag(V_p_mu))

    # Intervals
    z_critico <- stats::qnorm(0.975)
    ic_inferior <- estimativas_originais - z_critico * erros_padrao_orig
    ic_superior <- estimativas_originais + z_critico * erros_padrao_orig

    # Table setup
    TabelaR <- round(data.frame(
      MLE = as.numeric(estimativas_originais),
      SE = erros_padrao_orig,
      LI_05 = as.numeric(ic_inferior),
      LS_95 = as.numeric(ic_superior)
    ), 3)

    tabela_ic <- data.frame(Parameter = c('alpha', 'beta', 'a', 'b'), TabelaR)

    cat("\n--- Confidence Interval (95%) via Delta Method ---\n")
    print(tabela_ic)
  } else {
    cat("\n--- Confidence interval could not be calculated ---\n")
  }

  # 5. Quantile Residues
  # Note: Assumes that pwbs() is available in the package namespace
  theta <- stats::setNames(estimativas_originais, c('alpha', 'beta', 'a', 'b'))
  cdf_vals <- pwbs(y, theta['alpha'], theta['beta'], theta['a'], theta['b'])

  # Treatment for extreme values in CDF (avoids Inf in qnorm)
  cdf_vals <- pmax(pmin(cdf_vals, 1 - 1e-16), 1e-16)
  residuo <- stats::qnorm(cdf_vals)

  # Normality Test
  test_norm <- stats::shapiro.test(residuo)
  Wstats <- test_norm$statistic
  pw <- test_norm$p.value

  cat("\n---------- Diagnosis: Randomized Quantile Residuals ----------\n")
  cat("Shapiro-Wilk Test:\n")
  cat("Statistic:", round(Wstats, 3), "\n")
  cat("P-value:", round(pw, 3), "\n")

  if (pw <= 0.05) {
    cat("ATTENTION: Residuals are NOT normal (p < 0.05).\n")
  } else {
    cat("Residuals can be considered normal (p > 0.05).\n")
  }
  cat("--------------------------------------------------------------\n")

  # 6. Q-Q Chart
  plot_obj <- NULL
  if (plot) {
    plot_obj <- ggpubr::ggqqplot(residuo) +
      ggplot2::stat_qq(shape = 21, size = 2.5, col = 'blue2') +
      ggplot2::stat_qq_line(col = 'red', lwd = 0.7) +
      ggplot2::labs(
        title = "Q-Q Plot of model residuals",
        x = "Theoretical Quantiles",
        y = "Sample Residuals",
        caption = paste("Shapiro-Wilk: W =", round(Wstats, 3), ", p =", round(pw, 3))
      ) +
      ggplot2::theme_classic()

    print(plot_obj)
  }

  # Returns an invisible list for later use
  invisible(list(
    Estimates = theta,
    CI = tabela_ic,
    Metrics = tabela_info,
    Residuals = residuo,
    Normality_Test = test_norm,
    Plot = plot_obj
  ))
}
