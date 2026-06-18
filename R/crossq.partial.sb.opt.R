##' Returns critical values for the partial cross-quantilogram, based on the stationary bootstrap with the choice of the stationary-bootstrap parameter.
##'
##' This function generates critical values for for the partial cross-quantilogram,
##' using the stationary bootstrap in Politis and Romano (1994).
##' @title Stationary Bootstrap for the Partial Cross-Quantilogram dwith the choice of the stationary-bootstrap parameter
##' @param DATA The original data matrix
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k    A lag order
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level
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
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import np stats
##' @export

crossq.partial.sb.opt = function(DATA, vecA, k, Bsize, sigLev)
{
    ## optimal block length on all series (Politis-White), then delegate
    gamma = mean( 1/b.star(as.matrix(DATA))[,1] )
    crossq.partial.sb(DATA, vecA, k, gamma, Bsize, sigLev)

}  ## EoF
