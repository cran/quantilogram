% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Qstat.reg.sb.R
\name{Qstat.reg.sb}
\alias{Qstat.reg.sb}
\title{Stationary Bootstrap for Q statistics}
\usage{
Qstat.reg.sb(DATA1, DATA2, vecA, Psize, gamma, Bsize, sigLev)
}
\arguments{
\item{DATA1}{The original data set (1)}

\item{DATA2}{The original data set (2)}

\item{vecA}{A pair of two probabity values at which sample quantiles are estimated}

\item{Psize}{The maximum number of lags}

\item{gamma}{A parameter for the stationary bootstrap}

\item{Bsize}{The number of repetition of bootstrap}

\item{sigLev}{The statistical significance level}
}
\value{
The bootstrap critical values
}
\description{
Stationary Bootstrap procedure to generate critical values for both Box-Pierece and Ljung-Box type Q-statistics
}
\details{
This function returns critical values for for both Box-Pierece and Ljung-Box type Q-statistics through the statioanry bootstrap proposed by Politis and Romano (1994).
}
\examples{
data(sys.risk) 

## sample size
T = nrow(sys.risk)

## matrix for quantile regressions
## - 1st column: dependent variables
## - the rest:   regressors or predictors 
D1 = cbind(sys.risk[2:T,"Market"], sys.risk[1:(T-1),"Market"])
D2 = cbind(sys.risk[2:T,"JPM"], sys.risk[1:(T-1),"JPM"])

## probability levels
vecA = c(0.1, 0.2)

## setup for stationary bootstrap
gamma  = 1/10 ## bootstrap parameter depending on data
Bsize  = 5    ## small size, 5, for test 
sigLev = 0.05 ## significance level

## Q statistics with lags from 1 to 5, after quantile regression 
Qstat.reg.sb(D1, D2, vecA, 5, gamma, Bsize, sigLev)

}
\references{
Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
"The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
\emph{Journal of Econometrics}, 193(1), 251-270.

Politis, Dimitris N., and Joseph P. Romano. (1994).
"The stationary bootstrap."
\emph{Journal of the American Statistical Association} 89.428, pp.1303-1313.
}
\author{
Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
}
