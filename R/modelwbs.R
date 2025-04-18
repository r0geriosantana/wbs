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
#liklihood function
  LKH <- function(par) {
    alpha <- exp(par[1])
    beta <- exp(par[2])
    a <- exp(par[3])
    b <- exp(par[4])
    -sum(log(dwbs(y, alpha, beta, a, b)))
  }

  # Otimization
  resultado=nlminb(chute,LKH,lower = c(rep(-Inf,length(chute)),0), upper = c(rep(length(chute)),Inf))
  par_exp <- exp(resultado$par)
  theta <- setNames(par_exp, c('alpha', 'beta', 'a', 'b'))
  # AIC and BIC
  n <- length(y)
  n_param=length(chute)
  logL= resultado$objective #minus
  aic= 2 * logL + 2 * n_param
  bic= 2 * logL + n_param*log(n)


  cdf <- pwbs(y, theta['alpha'], theta['beta'], theta['a'], theta['b'])
  residuo <- qnorm(cdf)

  # Normality test (Shapiro-Wilk)
  test=shapiro.test(residuo)
  Wstats=test$statistic
  pw=test$p.value


  # Plot Q-Q plot if requested
  if(plot){
    plot_residual=ggqqplot(residuo) +
      stat_qq(shape = 21, size = 2.5, col = 'blue2') +
      stat_qq_line(col = 'red', lwd = 0.7) +
      labs(
        title = "Q-Q Plot of model residuals",
        x = "Theoretical Quantiles",
        y = "Sample Residuals",
        caption = paste("Shapiro-Wilk test, statistic =", round(Wstats, 4), ", p-value =", round(pw, 4))
      ) +
      theme_classic()
    print(plot_residual)
  }
  # List results
  result_list=list(
    MLE = round(theta, 3),
    AIC = round(aic, 3),
    BIC = round(bic, 3)
  )
result_list
}
