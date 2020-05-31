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
    ## size
    Tsize  = nrow(DATA)    ## =: T
    Nvar   = ncol(DATA)    ## =: #var
    Nsize  = Tsize - k

    ## =============================
    ## a pair of data points
    ## =============================
    ## original data, to draw {x_1t, x_2t-k}
    matD          = matrix(0, Nsize, Nvar)       ## N x 2
    matD[,1]      = DATA[(k+1):Tsize ,1     ,drop=FALSE]
    matD[,2:Nvar] = DATA[    1:Nsize ,2:Nvar,drop=FALSE]
    
    ##=========================================================================
    ## stationary bootstrap for cross-quantilogram
    ##=========================================================================
    ## container
    vecCRQ.B = matrix(0,Bsize,1) ## B x 1

    for (b in 1:Bsize){

        ##=============================
        ## Stationary Resampling
        ##=============================
        ## container 
        vecI = sb.index(Nsize, gamma)

        ## stationary resample
        matD.SB = matD[vecI,]       ## N x #var
        
        ## quantile hit
        matQhit.SB = q.hit(matD.SB, vecA)
        
        ## corss-quantilogram: not yet centered
        matHH.SB    = t(matQhit.SB) %*% matQhit.SB
        invHH.SB    = solve( matHH.SB)
        vecCRQ.B[b] = - invHH.SB[1,2] / sqrt( invHH.SB[1,1] * invHH.SB[2,2] )
        
    }

    ##=========================================================================
    ## Partial cross-quantilogram based on the original data
    ##=========================================================================
    RES     = crossq.partial(DATA, vecA, k)
    vParCRQ = RES$ParCRQ

    ## centering
    vecCRQ.cent = vecCRQ.B - vParCRQ  

    ## critical values
    vecCV    = matrix(0, 2, 1)
    vecCV[1] = quantile(vecCRQ.cent, (    sigLev / 2)) 
    vecCV[2] = quantile(vecCRQ.cent, (1 - sigLev / 2)) 

    ## results
    list(vecCV = vecCV, vParCRQ = vParCRQ)

}  ## EoF
