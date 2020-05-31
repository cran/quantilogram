##' The Data Set for Systemic Risk Analysis
##'
##' The data set contains the daily CRSP market value weighted index returns,
##' which are used as the market index returns in Brownless and Engle (2012), and
##' also includes daily stock returns on JP Morgan Chase (JPM), Goldman Sachs (GS)
##' and American International Group (AIG).
##' The sample period is from 2 Jan. 2001 to 30 Dec. 2011 with sample size 2,767.
##'
##' \itemize{
##'   \item date:    The time index (day)
##'   \item Market: The daily CRSP market value weighted incex returns
##'   \item JPM:    stock returns on JP Morgan Chase (JPM)
##'   \item GS:     stock returns on Goldman Sachs (GS)
##'   \item AIG:    stock returns on American International Group (AIG)
##' }
##'
##' @docType data
##' @keywords datasets
##' @name sys.risk
##' @usage data(sys.risk)
##' @format A data object with five variables
##' @references
##' Brownlees, Christian T., and Robert F. Engle. "Volatility, correlation and tails for systemic risk measurement." \emph{Available at SSRN} 1611229 (2012).
##'
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series."
##' \emph{Journal of Econometrics}, 193(1), 251-270.
##' 
##' 
"sys.risk"
