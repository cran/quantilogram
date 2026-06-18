##' Stationary Bootstrap procedure to generate critical values for both Box-Pierece and Ljung-Box type Q-statistics with the choice of the stationary-bootstrap parameter.
##'
##' This function returns critical values for for both Box-Pierece and Ljung-Box type
##' Q-statistics through the statioanry bootstrap proposed by Politis and Romano (1994).
##' To choose parameter for the statioanry bootstrap, this function
##' first obtaines the optimal value for each time serie using the
##' result provided by Politis and White (2004) and Patton, Politis and White (2004)
##' (The R-package, "np", written by  Hayfield and Racine is used).
##' Next, the average of the obtained values is used as the parameter value.
##' @title Stationary Bootstrap for Q statistics
##' @param DATA The original data
##' @param vecA A pair of two probabity values at which sample quantiles are estimated
##' @param Psize The maximum number of lags
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level
##' @return The bootstrap critical values
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series." \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' Patton, A., Politis, D. N., and White, H. (2009).
##' Correction to "Automatic block-length selection for the dependent bootstrap"
##' by D. Politis and H. White. \emph{Econometric Reviews}, 28(4), 372-375.
##'
##' Politis, D. N., and White, H. (2004).
##' "Automatic block-length selection for the dependent bootstrap."
##' \emph{Econometric Reviews}, 23(1), 53-70.
##'
##' Politis, Dimitris N., and Joseph P. Romano. (1994).
##' "The stationary bootstrap."
##' \emph{Journal of the American Statistical Association} 89.428: 1303-1313.
##'
##' @examples
##' data("sys.risk") ## data source
##' D = sys.risk[,c("Market", "JPM")] ## data: 2 variables
##'
##' # probability levels for the 2 variables
##' vecA = c(0.1, 0.5)
##'
##' ## setup for stationary bootstrap
##' Bsize  = 5    ## small size, 5, for test
##' sigLev = 0.05 ## significance level
##'
##' ## Q statistics with lags from 1 to5
##' Qstat.sb.opt(D, vecA, 5, Bsize, sigLev)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import np stats
##' @export

Qstat.sb.opt = function(DATA, vecA, Psize, Bsize, sigLev)
{
  ## optimal block length on the two series (Politis-White), then delegate
  gamma = mean( 1/b.star(as.matrix(DATA[, 1:2]))[,1] )
  Qstat.sb(DATA, vecA, Psize, gamma, Bsize, sigLev)

}  ## EOF
