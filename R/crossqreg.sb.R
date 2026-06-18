##' Returns critical values for the cross-quantilogram, based on the stationary bootstrap.
##'
##' This function generates critical values for for the cross-quantilogram,
##' using the stationary bootstrap in Politis and Romano (1994).
##' @title Stationary Bootstrap for the Cross-Quantilogram
##' @param DATA1 The original data matrix (T x p1)
##' @param DATA2 The original data matrix (T x p2)
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k    A lag order
##' @param gamma A parameter for the stationary bootstrap
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level
##' @return The boostrap critical values
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series." \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' Politis, Dimitris N., and Joseph P. Romano.
##' "The stationary bootstrap."
##' \emph{Journal of the American Statistical Association} 89.428 (1994): 1303-1313.
##'
##' @examples
##' data(sys.risk)
##'
##' ## sample size
##' T = nrow(sys.risk)
##'
##' ## matrix for quantile regressions
##' ## - 1st column: dependent variables
##' ## - the rest:   regressors or predictors
##' D1 = cbind(sys.risk[2:T,"Market"], sys.risk[1:(T-1),"Market"])
##' D2 = cbind(sys.risk[2:T,"JPM"], sys.risk[1:(T-1),"JPM"])
##'
##' ## probability levels
##' vecA = c(0.1, 0.2)
##'
##' ## setup for stationary bootstrap
##' gamma  = 1/10 ## bootstrap parameter depending on data
##' Bsize  = 5    ## small size, 5, for test
##' sigLev = 0.05 ## significance level
##'
##' ## cross-quantilogram with the lag of 5, after quantile regression
##' crossqreg.sb(D1, D2, vecA, 5, gamma, Bsize, sigLev)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats
##' @export

crossqreg.sb = function(DATA1, DATA2, vecA, k, gamma, Bsize, sigLev)
{
    ## size and the aligned pair {(y_1t, x_1t), (y_2,t-k, x_2,t-k)}
    Tsize = nrow(DATA1)
    Nsize = Tsize - k
    matD1 = DATA1[(k+1):Tsize, ]
    matD2 = DATA2[1:Nsize, ]                       ## already k-lagged

    ## one bootstrap cross-quantilogram (lag-0 cross-moment of the resampled hits)
    boot1 = function(b) {
        vecI  = sb.index(Nsize, gamma)
        matHH = crossprod(qreg.hit(matD1[vecI, ], matD2[vecI, ], vecA))
        matHH[1,2] / sqrt(matHH[1,1] * matHH[2,2])
    }
    vecCRQ.B = vapply(1:Bsize, boot1, numeric(1))

    ## cross-quantilogram on the original data, then centered critical values
    vCRQ  = crossqreg(DATA1, DATA2, vecA, k)
    vecCV = matrix(quantile(vecCRQ.B - vCRQ, c(sigLev/2, 1 - sigLev/2)), 2, 1)

    list(vecCV = vecCV, vCRQ = vCRQ)

}  ## EoF
