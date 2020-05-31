##' A function used to obtain partial cross-correlation function for a give lag order
##'
##' This function obtains the partial corss-correlation and the simple correlation.
##' To obtain the partial cross-correlation, this function uses the first column of
##' the input matrix and k-lagged values of the rest of the matrix.
##' @title  Partial Cross-correlation function
##' @param  matH A matrix with multiple columns (more than 3 columns)
##' @param  k    The lag order (integer)
##' @return Partial corss-correlation at k lags and the correlation statistics at k lags.
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang

corr.lag.partial = function(matH, k)
{
    ## size
    Tsize  = nrow(matH)    ## =: T
    Nvar   = ncol(matH)    ## =: #var
    Nsize  = Tsize - k

    ##  {H_1t, H_2t-k}
    matD          = matrix(0, Nsize, Nvar)
    matD[,1]      = as.matrix(matH[(k+1):Tsize,      1, drop=FALSE])
    matD[,2:Nvar] = as.matrix(matH[    1:Nsize, 2:Nvar, drop=FALSE]) 

    ## the following matrix contains inner-products of two vectors in H.
    matDD = t(matD) %*% matD

    ## cross-quantilogram of lag order k
    CRQ =  matDD[1,2] / sqrt( matDD[1,1] * matDD[2,2] ) ## 1 x 1

    ## partial quantilogram
    invDD  = solve(matDD)
    ParCRQ = - invDD[1,2] / sqrt( invDD[1,1] * invDD[2,2] )

    ## list
    list(CRQ = CRQ, ParCRQ = ParCRQ)

} ## EoF
