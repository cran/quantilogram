##' Returns critical values for the cross-quantilogram, based on the stationary bootstrap with the choice of the stationary-bootstrap parameter.
##'
##' This function generates critical values for for the cross-quantilogram,
##' using the stationary bootstrap in Politis and Romano (1994).
##' To choose parameter for the statioanry bootstrap, this function
##' first obtaines the optimal value for each time serie using the
##' result provided by Politis and White (2004) and Patton, Politis and White (2004)
##' (The R-package, "np", written by  Hayfield and Racine is used).
##' Next, the average of the obtained values is used as the parameter value.
##' @title Stationary Bootstrap for the Cross-Quantilogram with the choice of the stationary-bootstrap parameter
##' @param DATA The original data matrix
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k A lag order
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

crossq.sb.opt = function(DATA, vecA, k, Bsize, sigLev)
{
    ## size
    Tsize  = nrow(DATA)    ## =: T
    Nsize  = Tsize - k     ## =: N

    ## =============================
    ## a pair of data points
    ## =============================
    ## original data, to draw {x_1t, x_2t-k}
    matD     = matrix(0, Nsize, 2)              ## N x 2
    matD[,1] = as.matrix(DATA[(k+1):Tsize ,1, drop=FALSE]) ## N x 1
    matD[,2] = as.matrix(DATA[    1:Nsize ,2, drop=FALSE]) ## N x 1

    ##======================o
    ## optimal block size
    ##======================
    matB  = b.star(matD)                  ## use the sample being resampled
    gamma = mean( matB[,1, drop=FALSE] )  ## 1st colum = stationary boostrap

    ##=========================================================================
    ## stationary bootstrap for cross-quantilogram
    ##=========================================================================
    ## container
    vecCRQ.B = matrix(0,Bsize,1) ## B x 1

    for (b in 1:Bsize){

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
