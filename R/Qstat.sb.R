##' Stationary Bootstrap procedure to generate critical values for both Box-Pierece and Ljung-Box type Q-statistics
##'
##' This function returns critical values for for both Box-Pierece and Ljung-Box type Q-statistics through the statioanry bootstrap proposed by Politis and Romano (1994).
##' @title Stationary Bootstrap for Q statistics 
##' @param DATA The original data
##' @param vecA A pair of two probabity values at which sample quantiles are estimated 
##' @param Psize The maximum number of lags 
##' @param gamma A parameter for the stationary bootstrap
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level 
##' @return The bootstrap critical values 
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' Politis, Dimitris N., and Joseph P. Romano. (1994).
##' "The stationary bootstrap."
##' \emph{Journal of the American Statistical Association} 89.428, pp.1303-1313.
##'
##' @examples
##' data("sys.risk") ## data source 
##' D = sys.risk[,c("Market", "JPM")] ## data: 2 variables 
##'
##' # probability levels for the 2 variables 
##' vecA = c(0.1, 0.5)
##'
##' ## setup for stationary bootstrap
##' gamma  = 1/10 ## bootstrap parameter depending on data
##' Bsize  = 5    ## small size, 5, for test 
##' sigLev = 0.05 ## significance level
##'
##' ## Q statistics with lags from 1 to5
##' Qstat.sb(D, vecA, 5, gamma, Bsize, sigLev)
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats        
##' @export

Qstat.sb = function(DATA, vecA, Psize, gamma, Bsize, sigLev)
{
    ## size
    Tsize = nrow(DATA)    ## =: T    

    ##=========================================================================
    ## Q-stat based on the original data
    ##=========================================================================
    ## container
    vecQ.BP = matrix(0,Psize,1) ## P x 1
    vecQ.LB = matrix(0,Psize,1) ## P x 1

    ## cross-Q from 1 to P lags
    vecCRQ = crossq.max(DATA, vecA, Psize)  ## P x 1
    
    for (k in 1:Psize){
        
        ## Qstat
        RES        = Qstat(vecCRQ[1:k], Tsize)
        vecQ.BP[k] = RES$Q.BP
        vecQ.LB[k] = RES$Q.LB
    }

    ##=========================================================================
    ## stationary bootstrap for cross-quantilogram
    ##=========================================================================
    ## original data, to draw {x.1t, x.2t-1,...,x.2t-p}
    Nsize    = Tsize - Psize                  ## =: N
    matD     = matrix(0,Nsize,Psize+1)        ## N x p+1
    bigA     = matrix(0,(Psize+1), 1)         ## p+1 x 1
    matD[,1] = as.matrix(DATA[(Psize+1):Tsize,1,drop=F])   
    bigA[,1] = vecA[1]

    for (k in 1:Psize){
        matD[,(k+1)] = as.matrix(DATA[(Psize-k+1):(Tsize-k),2, drop=F])
        bigA[ (k+1)] = vecA[2]
    }


    ## containers
    matQ.BP = matrix(0,Bsize,Psize)  ## B x P: Box-Pierece Qstat
    matQ.LB = matrix(0,Bsize,Psize)  ## B x P: Ljung-Box Qstat

    for (b in 1:Bsize){
        
        ##=============================
        ## Stationary Resampling
        ##=============================
        ## selection index
        vecI = sb.index(Nsize, gamma)
        
        ## SB sample for {x.1t, x.2t-1, x.2t-2,..., x.2t-p}
        matD.SB = matD[vecI,]    

        ##=======================================
        ## Analysis based on stationary resample
        ##=======================================
        ##+++++++++++++++++++++++++++++
        ## Step 1: cross-quantilogram
        ##+++++++++++++++++++++++++++++
        vecCRQ.B = matrix(0,Psize, 1)

        ## quantile hit
        matQhit.SB = q.hit(matD.SB, bigA)
        
        for (k in 1:Psize){
            
            ## A pair of hits
            matH.SB     = cbind(matQhit.SB[,1], matQhit.SB[,(1+k)])
            matHH.SB    = t(matH.SB) %*% matH.SB
            vecCRQ.B[k] = matHH.SB[1,2] / sqrt( matHH.SB[1,1] * matHH.SB[2,2] ) ## 1x1
        }
        ##+++++++++++++++++++++++++++++
        ## Step 2: Q-stat 
        ##+++++++++++++++++++++++++++++
        ## centering
        vecTest = vecCRQ.B - vecCRQ
        
        ## Q stat
        for (k in 1:Psize){
            
            ## Q stat
            Res1         = Qstat(vecTest[1:k], Tsize)
            matQ.BP[b,k] = Res1$Q.BP
            matQ.LB[b,k] = Res1$Q.LB    
        }
    }

    ##=========================================================================
    ## critical values
    ##=========================================================================
    ## sort: column-wise
    vecCV.BP = matrix(0, Psize, 1)
    vecCV.LB = matrix(0, Psize, 1)

    ## critical values
    for (k in 1:Psize){

        vecCV.BP[k] = quantile(matQ.BP[,k], probs = (1 - sigLev) ) # 1 x 1
        vecCV.LB[k] = quantile(matQ.LB[,k], probs = (1 - sigLev) ) # 1 x 1
    }

    ## return
    list(vecQ.BP = vecQ.BP, vecCV.BP=vecCV.BP, vecQ.LB=vecQ.LB, vecCV.LB=vecCV.LB)

}  ## EOF

