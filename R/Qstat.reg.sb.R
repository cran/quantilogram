##' Stationary Bootstrap procedure to generate critical values for both Box-Pierece and Ljung-Box type Q-statistics
##'
##' This function returns critical values for for both Box-Pierece and Ljung-Box type Q-statistics through the statioanry bootstrap proposed by Politis and Romano (1994).
##' @title Stationary Bootstrap for Q statistics 
##' @param DATA1 The original data set (1)
##' @param DATA2 The original data set (2)
##' @param vecA A pair of two probabity values at which sample quantiles are estimated 
##' @param Psize The maximum number of lags 
##' @param gamma A parameter for the stationary bootstrap
##' @param Bsize The number of repetition of bootstrap
##' @param sigLev The statistical significance level 
##' @return The bootstrap critical values 
##' 
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
##' ## setup for stationary bootstrap
##' gamma  = 1/10 ## bootstrap parameter depending on data
##' Bsize  = 5    ## small size, 5, for test 
##' sigLev = 0.05 ## significance level
##' 
##' ## Q statistics with lags from 1 to 5, after quantile regression 
##' Qstat.reg.sb(D1, D2, vecA, 5, gamma, Bsize, sigLev)
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats        
##' @export

Qstat.reg.sb = function(DATA1, DATA2, vecA, Psize, gamma, Bsize, sigLev)
{
    ## size
    Tsize = nrow(DATA1)    ## =: T    
    Nsize = Tsize - Psize  ## =: N
    
    ##=========================================================================
    ## Q-stat based on the original data
    ##=========================================================================
    ## container
    vecQ.BP = matrix(0,Psize,1) ## P x 1
    vecQ.LB = matrix(0,Psize,1) ## P x 1

    ## cross-Q from 1 to P lags
    vecCRQ = crossqreg.max(DATA1, DATA2, vecA, Psize)  ## P x 1
    
    for (k in 1:Psize){
        
        ## Qstat
        RES        = Qstat(vecCRQ[1:k], Tsize)
        vecQ.BP[k] = RES$Q.BP
        vecQ.LB[k] = RES$Q.LB
    }

    ##=========================================================================
    ## stationary bootstrap for cross-quantilogram
    ##=========================================================================
    ## data 1 and 2
    matD1 = DATA1[(Psize+1):Tsize,] ## N x p1+1   
    
    ## containers
    matCRQ  = matrix(0,Bsize,Psize)  ## B x P: cross-quantilograms
    matQ.BP = matrix(0,Bsize,Psize)  ## B x P: Box-Pierece Qstat
    matQ.LB = matrix(0,Bsize,Psize)  ## B x P: Ljung-Box Qstat

    for (b in 1:Bsize){
        
        ##=============================
        ## Stationary Resampling
        ##=============================
        ## selection index
        vecI = sb.index(Nsize, gamma)
        ## - use the same index for all lags 
        
        ## SB sample from DATA1
        ## D1.SB: (y_11, x_11), ..., (y_1N,x_1N)
        matD1.SB = matD1[vecI,]    

        for (k in 1:Psize){
        
            ## SB sample from DATA2
            ## D2.SB: (y_2k, x_2k), ..., (y_1,N-k,x_1,N-k)
            matD2    = DATA2[(Psize-k+1):(Tsize-k),] ## N x p2+1   
            matD2.SB = matD2[vecI,]

            ##=======================================
            ## Analysis based on stationary resample
            ##=======================================
            matCRQ[b,k] = crossqreg(matD1.SB, matD2.SB, vecA, k)
        }

        ##+++++++++++++++++++++++++++++
        ## Step 2: Q-stat 
        ##+++++++++++++++++++++++++++++
        ## centering
        vecTest = matCRQ[b,] - vecCRQ ## P x 1
        
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
} ## EOF


