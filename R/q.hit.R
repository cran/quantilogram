##' Returns the matrix of quantil-hits
##'
##' This function generates the quantile hits given a vector of probabilty values.
##' The quantile hits are obtained for each column of an input matrix. 
##' @title Quantile Hit
##' @param DATA  A matrix that has time-series observations in its columns 
##' @param vecA  A vector of probabilty values at which sample quantiles are estimated 
##' @return      A matrix of quantile-hits
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @import stats 

q.hit = function(DATA, vecA)
{ 
  ## size
  Tsize  = nrow(DATA)    ## =: T
  Nvar   = ncol(DATA)    ## =: #var
  
  ## (1) sample quantile for each column
  vecQ = matrix(0,Nvar,1)           ## 2 x 1
  for (j in 1:Nvar){

      vecQ[j] = quantile(DATA[,j], probs = vecA[j])      
  }

  ## (2) Quantile Hit process with demean
  vecI    = matrix(1, Tsize, 1)
  matQhit = (DATA <= vecI %*% t(vecQ) ) - vecI %*% t(vecA)

  ## return
  return(matQhit)

}  ## EoF
