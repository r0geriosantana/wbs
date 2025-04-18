#' Kolmogorov-Smirnov Test for Weibull Birnbaum-Saunders Distribution
#'
#' @param y A numeric vector of data points.
#' @param alpha Shape parameter (alpha > 0) of the distribution.
#' @param beta Scale parameter (beta > 0) of the distribution.
#' @param a Shape parameter (a > 0) of the distribution.
#' @param b Shape parameter (b > 0) of the distribution.
#' @param alternative One of "two.sided" (default), "less", or "greater".
#' @param plot Logical; if TRUE, plots the CDFs.
#' @return A list with statistic and p-value.
#' @examples
#' alpha=1.5; beta=0.1; a=0.75;b=0.25
#' set.seed(152)
#' y <- rwbs(100,alpha,beta,a,b)
#' kswbs(y,alpha,beta,a,b,alternative = "two.sided", plot=FALSE)
#' @export
#'
kswbs <- function(y,alpha,beta,a,b, alternative = c("two.sided", "less", "greater"), plot = FALSE) {
  alternative <- match.arg(alternative)

  # Kolmogorov-Smirnov test
  ks_result <- ks.test(y, function(t) pwbs(t, alpha, beta, a, b), alternative = alternative)
  #statistics
  estatistica= ks_result$statistic
  #p-value
  pvalue = ks_result$p.value

    if(plot){
      empirical <- ecdf(y)
      x <- sort(y)
      empiric=empirical(x)
      theoretical <- sapply(x, function(t) pwbs(t, alpha, beta, a, b))
      plot=ggplot()+aes(x =x, y = theoretical) +
        geom_line(aes(x=x,y = empiric,color = "Empirical"), linetype = "dashed",size=0.7)+
        geom_line(aes(x=x,y = theoretical, color = "Theoretical"),size=0.7 )+
        scale_color_manual(
          name = "CDF",
          values = c("Empirical" = "blue", "Theoretical" = "red")
        ) +
        labs(title = "Kolmogorov-Smirnov Test",
             x = "t",
             y = "F(t)",
             caption = paste("Kolmogorov-Sminov test, D =",
                             round(estatistica, 4), ", p-value =", round(pvalue, 4)) ) +
        theme_classic()

      print(plot)  # Force the plot to display
    }

  list(statistic = estatistica, p.value =pvalue)
}
