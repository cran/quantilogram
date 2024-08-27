##' Plot of Cross-Quantilogram
##'
##' This function creates a plot of the cross-quantilogram with confidence intervals.
##' It computes the cross-quantilogram and its confidence intervals using stationary bootstrap,
##' then creates a ggplot visualization of the results.
##'
##' @param DATA A matrix of dimensions T x 2, where T is the number of observations.
##'             Column 1 contains the first variable and Column 2 contains the second variable.
##' @param vecA A numeric vector of quantiles for the first variable.
##' @param Kmax An integer representing the maximum lag to compute.
##' @param Bsize Bootstrap sample size for stationary bootstrap.
##' @param sigLev Significance level for confidence intervals. Default is 0.05 (95% confidence level).
##' @param vec.lag A vector of lag values (integer values). Not used in computation, only for plotting.
##' @param vec.CQ A numeric vector of cross-quantilogram values. Not used in computation, only for plotting.
##' @param mat.CI A matrix with two columns representing the lower and upper bounds of the confidence interval. Not used in computation, only for plotting.
##' @param y.min The minimum y-axis value. Default is -1.
##' @param y.max The maximum y-axis value. Default is 1.
##' @param ribbon_color Color for the confidence interval ribbon. Default is "gray".
##' @param ribbon_alpha Alpha (transparency) for the confidence interval ribbon. Default is 0.8.
##' @param bar_color Color for the quantilogram bars. Default is "black".
##' @param bar_width Width of the quantilogram bars. Default is 0.2.
##' @param title Plot title. Default is an empty string.
##' @param subtitle Plot subtitle. Default is NULL (no subtitle).
##'
#' @return A list containing two elements:
##'   \item{plot}{A ggplot object representing the cross-quantilogram plot over lags.}
##'   \item{df.res}{A data frame containing cross-quantilogram values and critical values. It includes the following columns:
##'     \itemize{
##'       \item lag: lag orders.
##'       \item crossQ: The cross-quantilogram values.
##'       \item CI_lower: The lower critical values for the confidence interval.
##'       \item CI_upper: The upper critical values for the confidence interval.
##'     }
##'   }
##'   
##' @return A list containing two elements:
##'   \item{plot}{A ggplot object representing the cross-quantilogram plot.}
##'   \item{df.res}{A data frame containing lag values, cross-quantilogram values, and confidence intervals.}
##'
##' @import ggplot2
##' @importFrom rlang .data
##'
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series." \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' @examples
##' \dontrun{
##' data("sys.risk")
##' DATA = sys.risk[,c("JPM", "Market")]
##' vecA = 0.05
##' Kmax = 20
##' Bsize = 200
##' result = crossq.plot(DATA, vecA, Kmax, Bsize)
##' print(result$plot)
##' }
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export
##'
crossq.plot = function(DATA, vecA, Kmax, Bsize, 
                       sigLev = 0.05, 
                       vec.lag, vec.CQ, mat.CI, 
                       y.min = -1, y.max = 1,
                       ribbon_color = "gray", ribbon_alpha = 0.8,
                       bar_color = "black", bar_width = 0.2,
                       title = "",
                       subtitle = NULL) {
  
  
  ## data frame: result container 
  df.res = data.frame(lag      = integer(Kmax), 
                      crossQ   = numeric(Kmax), 
                      CI_lower = numeric(Kmax), 
                      CI_upper = numeric(Kmax))
  
  ## cross-quantilogram + bootstrap 
  for (k in 1:Kmax) {
    RES = crossq.sb.opt(DATA, vecA, k, Bsize, sigLev)
    df.res$lag[k]      = k
    df.res$crossQ[k]   = RES$vCRQ
    df.res$CI_lower[k] = RES$vecCV[1]
    df.res$CI_upper[k] = RES$vecCV[2]
  }
  
  ## create x-axis label breaks 
  if (Kmax <= 5) {
    breaks = 1:Kmax
  } else {
    breaks = c(1, seq(5, Kmax, by = 5))
    if (tail(breaks, 1) != Kmax) breaks = c(breaks, Kmax)
  }
  
  ## plot over lags  
  p = ggplot(df.res, aes(x = .data$lag, y = .data$crossQ)) +
    geom_ribbon(aes(ymin = .data$CI_lower, 
                    ymax = .data$CI_upper), 
                fill = ribbon_color, alpha = ribbon_alpha) +
    geom_col(width = bar_width, fill = bar_color) +
    geom_hline(yintercept = 0) +
    labs(x = "Lag", y = "Quantilogram", 
         title = title, subtitle = subtitle) +
    ylim(y.min, y.max) +
    scale_x_continuous(breaks = breaks, labels = breaks) +
    theme_minimal()
  
  ## output 
  return(list(plot = p, df.res = df.res))
}
