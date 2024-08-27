##' Returns the cross-quantilogram 
##'
##' This function obtains the cross-quantilogram at the k lag order.
##' 
##' @title Cross-Quantilogram
##' @param DATA An input matrix of dimensions T x 2, where T is the number of observations.
##'             Column 1 contains the first variable and Column 2 contains the second variable.
##'             This function will apply a k-period lag to the second variable during computation.
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k    A lag order (integer)
##' @return Cross-Quantilogram
##'
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series." \emph{Journal of Econometrics}, 193(1), 251-270.
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
##' ## cross-quantilogram with the lag of 5
##' crossq.max(D, vecA, 5)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export

crossq = function(DATA, vecA, k)
{
    ## Quantile Hit process with demean
    matQhit = q.hit(DATA, vecA)

    ## correlation
    vCRQ = corr.lag(matQhit, k)

    ## cross-quantilogram
    return(vCRQ)  ## 1 x 1

}  ## EoF
