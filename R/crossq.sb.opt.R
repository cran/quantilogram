##' Returns critical values for the cross-quantilogram, based on the stationary bootstrap with the choice of the stationary-bootstrap parameter.
##'
##' This function generates critical values for for the cross-quantilogram,
##' using the stationary bootstrap in Politis and Romano (1994).
##' To choose parameter for the statioanry bootstrap, this function
##' first obtaines the optimal value for each time serie using the
##' result provided by Politis and White (2004) and Patton, Politis and White (2004)
##' (The R-package, "np", written by  Hayfield and Racine is used).
##' Next, the average of the obtained values is used as the parameter value.
##'
##' @title Stationary Bootstrap for the Cross-Quantilogram with the choice of the stationary-bootstrap parameter
##' @param DATA An input matrix of dimensions T x 2, where T is the number of observations.
##'             Column 1 contains the first variable and Column 2 contains the second variable.
##'             This function will apply a k-period lag to the second variable during computation.
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k A lag order
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level. Default is 0.05 (5% significance level).
##' @return The boostrap critical values
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
##' ## data source
##' data("sys.risk")
##'
##' ## data: 2 variables
##' D = sys.risk[,c("Market", "JPM")]
##'
##' # probability levels for the 2 variables
##' vecA = c(0.1, 0.5)
##'
##' ## setup for stationary bootstrap
##' Bsize  = 5    ## small size 5 for test
##' sigLev = 0.05 ## significance level
##'
##' ## cross-quantilogram with the lag of 5
##' crossq.sb.opt(D, vecA, 5, Bsize, sigLev)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import np stats
##' @export

crossq.sb.opt = function(DATA, vecA, k, Bsize, sigLev=0.05)
{
    ## optimal block length on the two series (Politis-White), then delegate
    gamma = mean( 1/b.star(as.matrix(DATA[,1:2]))[,1] )
    crossq.sb(DATA, vecA, k, gamma, Bsize, sigLev)

}  ## EoF
