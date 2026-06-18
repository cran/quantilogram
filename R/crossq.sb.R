##' Returns critical values for the cross-quantilogram, based on the stationary bootstrap.
##'
##' This function generates critical values for for the cross-quantilogram,
##' using the stationary bootstrap in Politis and Romano (1994).
##' @title Stationary Bootstrap for the Cross-Quantilogram
##' @param DATA An input matrix of dimensions T x 2, where T is the number of observations.
##'             Column 1 contains the first variable and Column 2 contains the second variable.
##'             This function will apply a k-period lag to the second variable during computation.
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k    A lag order
##' @param gamma A parameter for the stationary bootstrap
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level
##' @return The boostrap critical values
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' Politis, Dimitris N., and Joseph P. Romano. "The stationary bootstrap." \emph{Journal of the American Statistical Association} 89.428 (1994): 1303-1313.
##'
##' @examples
##' data("sys.risk") ## data source
##' D = sys.risk[,c("Market", "JPM")] ## data: 2 variables
##'
##' # probability levels for the 2 variables
##' vecA = c(0.1, 0.5)
##'
##' ## setup for stationary bootstrap
##' gamma  = 1/10 ## bootstrap parameter depending on data
##' Bsize  = 5    ## small size, 5, for test 
##' sigLev = 0.05 ## significance level
##'
##' ## cross-quantilogram with the lag of 5
##' crossq.sb(D, vecA, 5, gamma, Bsize, sigLev)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats
##' @export
crossq.sb = function(DATA, vecA, k, gamma, Bsize, sigLev)
{
    ## size and the aligned pair {x_1t, x_2,t-k}
    Tsize = nrow(DATA)
    Nsize = Tsize - k
    matD  = cbind(DATA[(k+1):Tsize, 1], DATA[1:Nsize, 2])   ## N x 2

    ## one bootstrap cross-quantilogram (lag-0 cross-moment of the resampled hits)
    boot1 = function(b) {
        matHH = crossprod(q.hit(matD[sb.index(Nsize, gamma), ], vecA))
        matHH[1,2] / sqrt(matHH[1,1] * matHH[2,2])
    }
    vecCRQ.B = vapply(1:Bsize, boot1, numeric(1))   ## B bootstrap values

    ## cross-quantilogram on the original data, then centered critical values
    vCRQ  = crossq(DATA, vecA, k)
    vecCV = matrix(quantile(vecCRQ.B - vCRQ, c(sigLev/2, 1 - sigLev/2)), 2, 1)

    list(vecCV = vecCV, vCRQ = vCRQ)

}  ## EoF
