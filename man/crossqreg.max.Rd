% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/crossqreg.max.R
\name{crossqreg.max}
\alias{crossqreg.max}
\title{Corss-Quantilogram up to a Given Lag Order}
\usage{
crossqreg.max(DATA1, DATA2, vecA, Kmax)
}
\arguments{
\item{DATA1}{An input matrix (T x p1)}

\item{DATA2}{An input matrix (T x p2)}

\item{vecA}{A pair of two probability values at which sample quantiles are estimated}

\item{Kmax}{The maximum lag order (integer)}
}
\value{
A vector of cross-quantilogram
}
\description{
The cross-quantilograms from 0 to a given lag order.
}
\details{
This function calculates the partial cross-quantilograms up to the lag order
users specify.
}
\references{
Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
"The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
\emph{Journal of Econometrics}, 193(1), 251-270.
}
\author{
Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
}
