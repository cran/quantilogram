##' A subfunction for the statioanry bootstrap
##'
##' This function resamples blocks of indicies with random block lengths.
##' This code follows the MATLAB file of the Oxford MFE Toolbox written by Kevin Sheppard.
##' @title Stationary Bootstrap Index 
##' @param Nsize The size of the stationary bootstrap resample 
##' @param gamma A parameter for the stationary boostrap.
##' @return A vector of indicies for the stationary bootstrap
##' @references
##' The Oxford MFE toolbox (http://www.kevinsheppard.com/wiki/MFE_Toolbox) by Kevin Sheppard
##' 
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang

sb.index = function(Nsize, gamma)
{
        
  ## container 
  vecI = matrix(0,Nsize,1) ## T x 1

  ## starting point
  vecI[1] = ceiling(Nsize * runif(n=1, min=0, max=1) ) 

  ## indicators: U[0,1] < gamma
  vecS = as.matrix(  ( runif(n=Nsize, min=0, max=1) < gamma) )   ## T x 1
  ## - if U[0,1] < gamma, then we pick up a new starting point.
  ## - if U[0,1] > gamma, then we stay. 

  vecU       = runif(n = sum(vecS), min=0, max=1)
  vecI[vecS] = ceiling( Nsize * vecU )

  for (i in 2:Nsize) {

    if ( vecS[i] == 0 ) {
      vecI[i] = vecI[i-1] + 1
    }

  }

  ## make it cicular: Exceeding
  vecE       = (vecI > Nsize)
  vecI[vecE] = vecI[vecE] - Nsize

  ## return
  return(vecI)

} ## EOF
