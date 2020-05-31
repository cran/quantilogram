##' Returns the partial cross-quantilogram 
##'
##' This function obtains the partial corss-quantilogram and the cross-quantilogram.
##' To obtain the partial cross-correlation given an input matrix, this function interacts
##' the values of the first column and the k-lagged values of the rest of the matrix.
##' @title Paritial Cross-Quantilogram
##' @param DATA1 An input matrix (T x p1)
##' @param DATA2 An input matrix (T x p2)
##' @param vecA A vector of probability values at which sample quantiles are estiamted 
##' @param k    The lag order 
##' @return The partial corss-quantilogram and the cross-quantilogram
##' 
##' @references 
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export 
crossqreg.partial = function(DATA1, DATA2, vecA, k)
{
    ## Important idea: make hit first
    ## Quantile Hit process with demean
    matH = qreg.hit(DATA1, DATA2, vecA)
    
    ## cross-quantilogram of lag order k
    RES = corr.lag.partial(matH, k)

    ## results 
    list(CRQ = RES$CRQ, ParCRQ = RES$ParCRQ  )

}  ## EoF
