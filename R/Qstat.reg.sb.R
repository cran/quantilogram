##' Stationary Bootstrap procedure to generate critical values for both Box-Pierece and Ljung-Box type Q-statistics
##'
##' This function returns critical values for for both Box-Pierece and Ljung-Box type Q-statistics through the statioanry bootstrap proposed by Politis and Romano (1994).
##' @title Stationary Bootstrap for Q statistics
##' @param DATA1 The original data set (1)
##' @param DATA2 The original data set (2)
##' @param vecA A pair of two probabity values at which sample quantiles are estimated
##' @param Psize The maximum number of lags
##' @param gamma A parameter for the stationary bootstrap
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level
##' @return The bootstrap critical values
##'
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
##' ## Q statistics with lags from 1 to 5, after quantile regression
##' Qstat.reg.sb(D1, D2, vecA, 5, gamma, Bsize, sigLev)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats
##' @export

Qstat.reg.sb = function(DATA1, DATA2, vecA, Psize, gamma, Bsize, sigLev)
{
    Tsize = nrow(DATA1)
    Nsize = Tsize - Psize

    ## Box-Pierce / Ljung-Box Q-statistics on the original data, lags 1..Psize
    vecCRQ  = crossqreg.max(DATA1, DATA2, vecA, Psize)
    Q       = vapply(1:Psize, function(k) unlist(Qstat(vecCRQ[1:k], Tsize)), numeric(2))
    vecQ.BP = matrix(Q[1, ], Psize, 1)
    vecQ.LB = matrix(Q[2, ], Psize, 1)

    ## DATA1 aligned at t; DATA2 pre-lagged by each k (one resample index for all lags)
    matD1 = DATA1[(Psize+1):Tsize, ]
    matD2 = lapply(1:Psize, function(k) DATA2[(Psize-k+1):(Tsize-k), ])

    ## one bootstrap draw -> centered Box-Pierce & Ljung-Box Q at lags 1..Psize
    boot1 = function(b) {
        vecI = sb.index(Nsize, gamma)
        d1   = matD1[vecI, ]
        cq   = vapply(1:Psize, function(k) crossqreg(d1, matD2[[k]][vecI, ], vecA, 0), numeric(1))
        vapply(1:Psize, function(k) unlist(Qstat((cq - vecCRQ)[1:k], Tsize)), numeric(2))
    }
    B = vapply(1:Bsize, boot1, matrix(0, 2, Psize))   ## 2 x Psize x Bsize

    cv = function(i) matrix(apply(B[i, , ], 1, quantile, probs = 1 - sigLev), Psize, 1)
    list(vecQ.BP = vecQ.BP, vecCV.BP = cv(1), vecQ.LB = vecQ.LB, vecCV.LB = cv(2))

} ## EOF
