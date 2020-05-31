##' Returns the cross-quantilogram 
##'
##' This function obtains the cross-quantilogram at the k lag order.
##' @title Cross-Quantilogram
##' @param DATA1 An input matrix (T x p1)
##' @param DATA2 An input matrix (T x p2)
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k    A lag order (integer)
##' @return Cross-Quantilogram
##'
##' @references 
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' Koenker, R., and Bassett Jr, G. (1978). 
##' "Regression quantiles." \emph{Econometrica}, 46(1), 33-50.
##'
##' @examples
##' ## data source 
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
##' ## cross-quantilogram with the lag of 5, after quantile regression 
##' crossqreg(D1, D2, vecA, 5)
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export

crossqreg = function(DATA1, DATA2, vecA, k)
{
    ## Quantile Hit process with demean
    matQhit = qreg.hit(DATA1, DATA2, vecA)

    ## correlation
    vCRQ = corr.lag(matQhit, k)

    ## cross quantilogram 
    return(vCRQ)  ## 1 x 1

}  ## EoF
