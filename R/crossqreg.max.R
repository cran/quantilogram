##' The cross-quantilograms from 0 to a given lag order.
##'
##' This function calculates the partial cross-quantilograms up to the lag order
##' users specify.
##' @title Corss-Quantilogram up to a Given Lag Order
##' @param DATA1 An input matrix (T x p1)
##' @param DATA2 An input matrix (T x p2)
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param Kmax The maximum lag order (integer)
##' @return A vector of cross-quantilogram
##' 
##' @references 
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export 
crossqreg.max = function(DATA1, DATA2, vecA, Kmax)
{
    ## Important idea: make hit first
    ## Quantile Hit process with demean
    matQhit = qreg.hit(DATA1, DATA2, vecA)
  
    ## for each lag
    vecCRQ = matrix(0, (Kmax+1), 1)  ## K+1 x 1
    
    for (k in 1:(Kmax+1)){
        ## cross-quantilogram of lag order k
        vecCRQ[k] =  corr.lag(matQhit, k)
    }

    ## vector of cross-quantilograms
    return(vecCRQ)  ## K+1 x 1

}  ## EoF
