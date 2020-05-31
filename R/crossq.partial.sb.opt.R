##' Returns critical values for the partial cross-quantilogram, based on the stationary bootstrap with the choice of the stationary-bootstrap parameter. 
##'
##' This function generates critical values for for the partial cross-quantilogram,
##' using the stationary bootstrap in Politis and Romano (1994).
##' @title Stationary Bootstrap for the Partial Cross-Quantilogram dwith the choice of the stationary-bootstrap parameter
##' @param DATA The original data matrix 
##' @param vecA A pair of two probability values at which sample quantiles are estimated 
##' @param k    A lag order 
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level 
##' @return The boostrap critical values
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series." \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' Patton, A., Politis, D. N., and White, H. (2009).
##' Correction to "Automatic block-length selection for the dependent bootstrap"
##' by D. Politis and H. White. \emph{Econometric Reviews}, 28(4), 372-375.
##'
##' Politis, D. N., and White, H. (2004).
##' "Automatic block-length selection for the dependent bootstrap."
##' \emph{Econometric Reviews}, 23(1), 53-70.
##'
##' Politis, Dimitris N., and Joseph P. Romano. (1994).
##' "The stationary bootstrap."
##' \emph{Journal of the American Statistical Association} 89.428: 1303-1313.
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import np stats
##' @export 

crossq.partial.sb.opt = function(DATA, vecA, k, Bsize, sigLev)
{
    ## size
    Tsize  = nrow(DATA)    ## =: T
    Nvar   = ncol(DATA)    ## =: #var
    Nsize  = Tsize - k     ## =: N

    ## =============================
    ## a pair of data points
    ## =============================
    ## original data, to draw {x_1t, x_2t-k}
    matD          = matrix(0, Nsize, Nvar)                ## N x 2
    matD[,1]      = DATA[(k+1):Tsize ,1     ,drop=FALSE]
    matD[,2:Nvar] = DATA[    1:Nsize ,2:Nvar,drop=FALSE]

    ##======================o
    ## optimal block size
    ##======================
    matB  = b.star(matD)                  ## use the sample being resampled
    gamma = mean( matB[,1, drop=FALSE] )  ## 1st colum = stationary boostrap
    
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
