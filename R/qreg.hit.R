##' Returns the matrix of quantil-hits
##'
##' This function generates the quantile hits based on quantile regression,
##' given a vector of probabilty values. The quantile regressions are esimated
##' for each matrix of data and a pair of quantile hits are produced.
##' @title Quantile Hit
##' @param DATA1 An input matrix (T x p1+1) with the first column of the dependent varaible and the the rest of columns with regressors
##' @param DATA2 An input matrix (T x p2+1) with the first column of the dependent varaible and the the rest of columns with regressors
##' @param vecA  A vector of probabilty values at which sample quantiles are estimated 
##' @return      A matrix of quantile-hits
##' @references 
##' Koenker, R., and Bassett Jr, G. (1978). 
##' "Regression quantiles." Econometrica, 46(1), 33-50.
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import quantreg

qreg.hit = function(DATA1, DATA2, vecA)
{
    ## size
    Tsize  = nrow(DATA1)    ## =: T
    p1     = ncol(DATA1)
    p2     = ncol(DATA2)
    
    ## Quantile Regression
    qfit1 = rq(DATA1[,1] ~ DATA1[,2:p1], vecA[1])
    qfit2 = rq(DATA2[,1] ~ DATA2[,2:p2], vecA[2])

    ## Residuals
    vecRes1 = qfit1$residuals
    vecRes2 = qfit2$residuals
    matRes  = data.matrix( cbind(vecRes1, vecRes2) )

    ##  Quantile Hit process with demean
    vecI    = matrix(1, Tsize, 1)        ## T x 1
    mat0    = matrix(0, nrow = Tsize, 2) ## T x 2   
    matQhit = (matRes <= mat0 ) - vecI %*% t(vecA)

    ## return
    return(matQhit)
    
}  ## EoF
