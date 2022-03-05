##' The partial cross-quantilograms from 1 to a given lag order.
##'
##' This function calculates the partial cross-quantilograms up to the lag order
##' users specify.
##' @title Partial Corss-Quantilogram upto a given lag order
##' @param DATA An input matrix 
##' @param vecA A vector of probability values at which sample quantiles are estimated
##' @param Kmax The maximum lag order (integer)
##' @return A vector of cross-quantilogram and a vector of partial cross-quantilograms
##'
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series." \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' @examples
##' ## data source 
##' data("sys.risk") 
##'
##' ## data with 3 variables 
##' D = sys.risk[,c("Market", "JPM", "VIX")]
##'
##' ## probablity levels for the 3 variables 
##' vecA = c(0.1, 0.1, 0.1)
##'
##' ## partial cross-quantilogram with lags from 1 to 5
##' crossq.max.partial(D, vecA, 5)
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export 
crossq.max.partial = function(DATA, vecA, Kmax)
{
    ## Quantile Hit process with demean
    matH = q.hit(DATA, vecA)
    
    ## (3) for each lag
    vecCRQ    = matrix(0, Kmax, 1)
    vecParCRQ = matrix(0, Kmax, 1)

    for (k in 1:Kmax){

        ## cross-quantilogram of lag order k
        RES          = corr.lag.partial(matH, k)
        vecCRQ[k]    = RES$CRQ
        vecParCRQ[k] = RES$ParCRQ

    }

    ## return 
    list(CRQ = vecCRQ, ParCRQ = vecParCRQ)

}  ## EoF
