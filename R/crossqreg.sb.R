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
##' Bsize  = 5    ## small size 10 for test 
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
    ## size
    Tsize  = nrow(DATA1)    ## =: T
    Nsize  = Tsize - k      ## =: N

    ## =============================
    ## a pair of data points
    ## =============================
    ## original data, to draw {(y_1t, x_1t) (y_2t-k, x_2t-k)}
    matD1 = DATA1[(k+1):Tsize,]  ## N x k2
    matD2 = DATA2[   1 :Nsize,]  ## N x k1 <== This is k-lagged value 
    
    ##=========================================================================
    ## stationary bootstrap for cross-quantilogram
    ##=========================================================================
    ## container
    vecCRQ.B = matrix(0,Bsize,1) ## B x 1

    for (b in 1:Bsize){

        ##=============================
        ## Stationary Resampling
        ##=============================
        vecI  = sb.index(Nsize, gamma) ## resample index 
        matD1.SB = matD1[vecI,]
        matD2.SB = matD1[vecI,]

        ## quantile hit
        matQhit.SB = qreg.hit(matD1.SB, matD2.SB, vecA)
        
        ## corss-quantilogram: not yet centered
        matHH.SB    = t(matQhit.SB) %*% matQhit.SB
        vecCRQ.B[b] = matHH.SB[1,2] / sqrt( matHH.SB[1,1] * matHH.SB[2,2] ) ## 1 x 1
    }

    ##=========================================================================
    ## cross-quantilogram based on the original data
    ##=========================================================================
    ## corss-quantilogram
    vCRQ = crossqreg(DATA1, DATA2, vecA, k)

    ## centering
    vecCRQ.cent = vecCRQ.B - vCRQ  

    ## critical values
    vecCV    = matrix(0, 2, 1)
    vecCV[1] = quantile(vecCRQ.cent, (    sigLev / 2)) 
    vecCV[2] = quantile(vecCRQ.cent, (1 - sigLev / 2)) 

    ## results
    list(vecCV = vecCV, vCRQ = vCRQ)

}  ## EoF
