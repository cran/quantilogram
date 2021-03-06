% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Qstat.R
\name{Qstat}
\alias{Qstat}
\title{Q-statistics}
\usage{
Qstat(vecTest, Tsize)
}
\arguments{
\item{vecTest}{A vector of test statistics ordered with respect the number of lags}

\item{Tsize}{A original sample size}
}
\value{
the Box-Pierece and Ljung-Box statistics
}
\description{
Te Box-Pierece and Ljung-Box type Q-statistics
}
\details{
This function returns Box-Pierece and Ljung-Box type Q-statistics
}
\references{
Box, G. EP, and D. A. Pierce. (1970).
"Distribution of residual autocorrelations in autoregressive-integrated moving average time series models."
\emph{Journal of the American Statistical Association} 65.332, pp.1509-1526.

Ljung, G. M., and G. EP Box. (1978).
"On a measure of lack of fit in time series models."
\emph{Biometrika} 65.2, pp.297-303.
}
\author{
Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
}
