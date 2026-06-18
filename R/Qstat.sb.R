##' Stationary Bootstrap procedure to generate critical values for both Box-Pierece and Ljung-Box type Q-statistics
##'
##' This function returns critical values for for both Box-Pierece and Ljung-Box type Q-statistics through the statioanry bootstrap proposed by Politis and Romano (1994).
##' @title Stationary Bootstrap for Q statistics
##' @param DATA The original data
##' @param vecA A pair of two probabity values at which sample quantiles are estimated
##' @param Psize The maximum number of lags
##' @param gamma A parameter for the stationary bootstrap
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level
##' @return The bootstrap critical values
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' Politis, Dimitris N., and Joseph P. Romano. (1994).
##' "The stationary bootstrap."
##' \emph{Journal of the American Statistical Association} 89.428, pp.1303-1313.
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
##' ## Q statistics with lags from 1 to5
##' Qstat.sb(D, vecA, 5, gamma, Bsize, sigLev)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats
##' @export

Qstat.sb = function(DATA, vecA, Psize, gamma, Bsize, sigLev)
{
    Tsize = nrow(DATA)
    Nsize = Tsize - Psize

    ## Box-Pierce / Ljung-Box Q-statistics on the original data, lags 1..Psize
    vecCRQ  = crossq.max(DATA, vecA, Psize)
    Q       = vapply(1:Psize, function(k) unlist(Qstat(vecCRQ[1:k], Tsize)), numeric(2))
    vecQ.BP = matrix(Q[1, ], Psize, 1)
    vecQ.LB = matrix(Q[2, ], Psize, 1)

    ## omnibus tuple {x_1t, x_2,t-1, ..., x_2,t-Psize} and its quantile levels
    matD = cbind(DATA[(Psize+1):Tsize, 1],
                 sapply(1:Psize, function(k) DATA[(Psize-k+1):(Tsize-k), 2]))
    bigA = c(vecA[1], rep(vecA[2], Psize))

    ## one bootstrap draw -> centered Box-Pierce & Ljung-Box Q at lags 1..Psize
    boot1 = function(b) {
        H  = q.hit(matD[sb.index(Nsize, gamma), ], bigA)
        h0 = H[, 1]; hk = H[, -1, drop = FALSE]
        cq = as.vector(crossprod(h0, hk)) / (sqrt(sum(h0^2)) * sqrt(colSums(hk^2)))
        vapply(1:Psize, function(k) unlist(Qstat((cq - vecCRQ)[1:k], Tsize)), numeric(2))
    }
    B = vapply(1:Bsize, boot1, matrix(0, 2, Psize))   ## 2 x Psize x Bsize

    cv = function(i) matrix(apply(B[i, , ], 1, quantile, probs = 1 - sigLev), Psize, 1)
    list(vecQ.BP = vecQ.BP, vecCV.BP = cv(1), vecQ.LB = vecQ.LB, vecCV.LB = cv(2))

}  ## EOF
