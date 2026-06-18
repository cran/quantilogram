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
##' @param y.min The minimum y-axis value. Default is -1.
##' @param y.max The maximum y-axis value. Default is 1.
##' @param ribbon_color Color for the confidence interval ribbon. Default is "gray".
##' @param ribbon_alpha Alpha (transparency) for the confidence interval ribbon. Default is 0.8.
##' @param bar_color Color for the quantilogram bars. Default is "black".
##' @param bar_width Width of the quantilogram bars. Default is 0.2.
##' @param title Plot title. Default is an empty string.
##' @param subtitle Plot subtitle. Default is NULL (no subtitle).
##' @param gamma Stationary-bootstrap parameter (the mean block length is 1/gamma). If NULL (default),
##'              it is chosen once by the Politis-White rule via \code{np::b.star} on the two series.
##'              Supplying a value avoids the dependency on the \code{np} package.
##'
##' @return A list containing two elements:
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
##' @import ggplot2
##' @importFrom rlang .data
##' @importFrom stats quantile
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
                       y.min = -1, y.max = 1,
                       ribbon_color = "gray", ribbon_alpha = 0.8,
                       bar_color = "black", bar_width = 0.2,
                       title = "",
                       subtitle = NULL,
                       gamma = NULL) {

  ## allow a single common quantile for both series
  if (length(vecA) == 1) {
    vecA = rep(vecA, 2)
  }

  X  = as.matrix(DATA[, 1:2])
  Tn = nrow(X)
  if (Kmax < 1 || Kmax >= Tn) {
    stop("Kmax must satisfy 1 <= Kmax < nrow(DATA)")
  }

  ## point estimates for lags 1..Kmax (full-sample quantiles, as in crossq)
  vecCRQ = crossq.max(DATA, vecA, Kmax)        ## Kmax x 1

  ## block-length parameter, chosen once
  if (is.null(gamma)) {
    gamma = mean(1 / np::b.star(X)[, 1])
  }

  ## All lags share one resample per draw: build the (y1_t, y2_{t-1}, ..., y2_{t-Kmax})
  ## tuple once, then read every lag off the resampled hits (the omnibus structure).
  Nn   = Tn - Kmax
  matD = matrix(0, Nn, Kmax + 1)
  matD[, 1] = X[(Kmax + 1):Tn, 1]
  for (k in 1:Kmax) {
    matD[, k + 1] = X[(Kmax - k + 1):(Tn - k), 2]
  }
  bigA = c(vecA[1], rep(vecA[2], Kmax))

  ## bootstrap: resample once per draw, all lags from one hit matrix
  boot = matrix(0, Bsize, Kmax)
  for (b in 1:Bsize) {
    vecI = sb.index(Nn, gamma)
    H    = q.hit(matD[vecI, ], bigA)            ## Nn x (Kmax+1)
    h0   = H[, 1]
    hk   = H[, -1, drop = FALSE]
    boot[b, ] = as.vector(crossprod(h0, hk)) / (sqrt(sum(h0^2)) * sqrt(colSums(hk^2)))
  }

  ## centred critical values per lag
  df.res = data.frame(lag      = 1:Kmax,
                      crossQ   = as.vector(vecCRQ),
                      CI_lower = numeric(Kmax),
                      CI_upper = numeric(Kmax))
  for (k in 1:Kmax) {
    cv = quantile(boot[, k] - vecCRQ[k], probs = c(sigLev / 2, 1 - sigLev / 2), na.rm = TRUE)
    df.res$CI_lower[k] = cv[1]
    df.res$CI_upper[k] = cv[2]
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
