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
    ## size
    Tsize  = nrow(DATA)    ## =: T
    Nsize  = Tsize - k     ## =: N

    ## =============================
    ## a pair of data points
    ## =============================
    ## original data, to draw {x_1t, x_2t-k}
    matD     = matrix(0, Nsize, 2)              ## N x 2
    matD[,1] = as.matrix(DATA[(k+1):Tsize, 1, drop=FALSE]) ## N x 1
    matD[,2] = as.matrix(DATA[    1:Nsize, 2, drop=FALSE]) ## N x 1

    ##=========================================================================
    ## stationary bootstrap for cross-quantilogram
    ##=========================================================================
    ## container
    vecCRQ.B = matrix(0,Bsize,1) ## B x 1

    for (b in 1:Bsize){

        ##=============================
        ## Stationary Resampling
        ##=============================
        ## container
        vecI = sb.index(Nsize, gamma)

        ## stationary resample
        matD.SB = matD[vecI,]       ## N x 2

        ## quantile hit
        matQhit.SB = q.hit(matD.SB, vecA)

        ## corss-quantilogram: not yet centered
        matHH.SB    = t(matQhit.SB) %*% matQhit.SB
        vecCRQ.B[b] = matHH.SB[1,2] / sqrt( matHH.SB[1,1] * matHH.SB[2,2] ) ## 1 x 1
    }

    ##=========================================================================
    ## cross-quantilogram based on the original data
    ##=========================================================================
    ## corss-quantilogram
    vCRQ = crossq(DATA, vecA, k)

    ## centering
    vecCRQ.cent = vecCRQ.B - vCRQ

    ## critical values
    vecCV    = matrix(0, 2, 1)
    vecCV[1] = quantile(vecCRQ.cent, (    sigLev / 2))
    vecCV[2] = quantile(vecCRQ.cent, (1 - sigLev / 2))

    ## results
    list(vecCV = vecCV, vCRQ = vCRQ)

}  ## EoF
