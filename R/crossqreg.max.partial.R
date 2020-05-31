##' The partial cross-quantilograms from 1 to a given lag order.
##'
##' This function calculates the partial cross-quantilograms up to the lag order
##' users specify.
##' @title Partial Corss-Quantilogram upto a given lag order
##' @param DATA1 An input matrix (T x p1)
##' @param DATA2 An input matrix (T x p2)
##' @param vecA A vector of probability values at which sample quantiles are estimated
##' @param Kmax The maximum lag order (integer)
##' @return A vector of cross-quantilogram and a vector of partial cross-quantilograms
##' 
##' @references 
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export
 
crossqreg.max.partial = function(DATA1, DATA2, vecA, Kmax)
{
    ## Quantile Hit process with demean
    matH = qreg.hit(DATA1, DATA2, vecA)
    
    ## (3) for each lag
    vecCRQ    = matrix(0, Kmax, 1)
    vecParCRQ = matrix(0, Kmax, 1)

    for (k in 1:Kmax){
        ## cross-quantilogram of lag order k
        RES          = corr.lag.partial(matH, k)
        vecCRQ[k]    = RES$CRQ
        vecParCRQ[k] = RES$ParCRQ
    }

    ## list
    list(CRQ = vecCRQ, ParCRQ = vecParCRQ)

}  ## EoF
