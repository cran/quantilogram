##' Te Box-Pierece and Ljung-Box type Q-statistics 
##'
##' This function returns Box-Pierece and Ljung-Box type Q-statistics 
##' @title Q-statistics 
##' @param vecTest A vector of test statistics ordered with respect the number of lags
##' @param Tsize   A original sample size
##' @return the Box-Pierece and Ljung-Box statistics 
##' @references
##' Box, G. EP, and D. A. Pierce. (1970).
##' "Distribution of residual autocorrelations in autoregressive-integrated moving average time series models."
##' \emph{Journal of the American Statistical Association} 65.332, pp.1509-1526.
##'
##' Ljung, G. M., and G. EP Box. (1978).
##' "On a measure of lack of fit in time series models."
##' \emph{Biometrika} 65.2, pp.297-303.
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang

Qstat = function(vecTest, Tsize)
{
    ## sizes 
    Psize = length(vecTest)
    vecP  = matrix( seq((Tsize - 1),(Tsize-Psize)), Psize, 1)
    invP  = vecP ^ (-1)
    vecT2 = matrix( (vecTest^2), Psize,1)

    ##================================
    ## Box-Pierece test
    ##================================
    Q.BP = Tsize * sum( vecT2 ) 

    ##=================================
    ## Ljung-Box type Q statitics
    ##=================================
    Q.LB  = Tsize * (Tsize + 2) * sum( (vecT2 * invP) )

    ## return
    list(Q.BP = Q.BP, Q.LB = Q.LB) 

}  ## END of function

