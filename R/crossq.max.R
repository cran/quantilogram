##' The cross-quantilograms from 1 to a given lag order.
##'
##' This function calculates the partial cross-quantilograms up to the lag order
##' users specify.
##' @title Corss-Quantilogram up to a Given Lag Order
##' @param DATA An input matrix 
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param Kmax The maximum lag order (integer)
##' @return A vector of cross-quantilogram
##'
##' @references 
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
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
##' ## cross-quantilogram with lags between 1 and 5
##' crossq.max(D, vecA, 5)
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export 
crossq.max = function(DATA, vecA, Kmax)
{
    ## Important idea: make hit first
    ## Quantile Hit process with demean
    matQhit = q.hit(DATA, vecA)
    
    ## for each lag
    vecCRQ = matrix(0, Kmax, 1)  ## K x 1
    
    for (k in 1:Kmax){
        
        ## cross-quantilogram of lag order k
        vecCRQ[k] =  corr.lag(matQhit, k)
    }

    ## 
    return(vecCRQ)  ## K x 1

}  ## EoF
