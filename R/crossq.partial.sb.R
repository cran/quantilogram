##' Returns critical values for the partial cross-quantilogram, based on the stationary bootstrap.
##'
##' This function generates critical values for for the partial cross-quantilogram,
##' using the stationary bootstrap in Politis and Romano (1994).
##' @title Stationary Bootstrap for the Partial Cross-Quantilogram
##' @param DATA The original data matrix
##' @param vecA A pair of two probability values at which sample quantiles are estimated
##' @param k    A lag order
##' @param gamma A parameter for the stationary bootstrap
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level
##' @return The boostrap critical values
##' @references
##' Politis, Dimitris N., and Joseph P. Romano. "The stationary bootstrap." \emph{Journal of the American Statistical Association} 89.428 (1994): 1303-1313.
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats
##' @export

crossq.partial.sb = function(DATA, vecA, k, gamma, Bsize, sigLev)
{
    ## size and the aligned data (col 1 at t, cols 2..Nvar at t-k)
    Tsize = nrow(DATA)
    Nvar  = ncol(DATA)
    Nsize = Tsize - k
    matD  = cbind(DATA[(k+1):Tsize, 1], as.matrix(DATA[1:Nsize, 2:Nvar]))   ## N x Nvar

    ## one bootstrap partial cross-quantilogram (NA if the moment matrix is singular)
    boot1 = function(b) {
        matHH = crossprod(q.hit(matD[sb.index(Nsize, gamma), ], vecA))
        invHH = tryCatch(solve(matHH), error = function(e) matrix(NA_real_, Nvar, Nvar))
        -invHH[1,2] / sqrt(invHH[1,1] * invHH[2,2])
    }
    vecCRQ.B = vapply(1:Bsize, boot1, numeric(1))

    ## partial cross-quantilogram on the original data, then centered critical values
    vParCRQ = crossq.partial(DATA, vecA, k)$ParCRQ
    vecCV   = matrix(quantile(vecCRQ.B - vParCRQ, c(sigLev/2, 1 - sigLev/2), na.rm = TRUE), 2, 1)

    list(vecCV = vecCV, vParCRQ = vParCRQ)

}  ## EoF
