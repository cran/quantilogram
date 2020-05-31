##' The correlation statistics for a given lag order 
##'
##' The function obtains the simple correlation statistics. The values in the first column of input matrix
##' is interacted with the k-lagged values in the second column.
##' @title Correlation Function
##' @param matH The matrix with the column size of 2 
##' @param k    The lag order (integer)
##' @return Correlation
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang

corr.lag = function(matH, k)
{
    ## size
    Tsize = nrow(matH)    ## =: T
    Nsize = Tsize - k     ## =: N

    ##  {H_1t, H_2t-k}
    matD     = matrix(0, Nsize, 2)  ## N x 2
    matD[,1] = as.matrix(matH[(k+1):Tsize, 1])  
    matD[,2] = as.matrix(matH[    1:Nsize, 2])
    
    ## the following matrix contains inner-products of two vectors in H.
    matDD  = t(matD) %*% matD   ## 2 x 2

    ## cross-quantilogram of lag order k
    CRQ =  matDD[1,2] / sqrt( matDD[1,1] * matDD[2,2] ) ## 1 x 1
    
    ## list
    return(CRQ)

} ## EoF
