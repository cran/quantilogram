#' Quantilogram Analysis Tools
#'
#' @description
#' This package provides a comprehensive set of tools for quantilogram analysis in R.
#' It includes functions for computing and visualizing cross-quantilograms, which are
#' useful for analyzing dependence structures in financial time series data.
#' The package implements methods described in Han et al. (2016) for measuring
#' quantile dependence and testing directional predictability between time series.
#'
#' @details
#' The package's functions can be categorized into several groups:
#'
#' \strong{Core Quantilogram Functions:}
#' \itemize{
#'   \item \code{\link{crossq}}: Compute basic cross-quantilogram
#'   \item \code{\link{crossq.sb}}: Cross-quantilogram with stationary bootstrap
#'   \item \code{\link{crossq.sb.opt}}: Optimized cross-quantilogram with bootstrap
#' }
#'
#' \strong{Visualization Functions:}
#' \itemize{
#'   \item \code{\link{crossq.heatmap}}: Create heatmap visualization of cross-quantilograms
#'   \item \code{\link{crossq.plot}}: Plot method for crossq objects
#' }
#'
#' \strong{Advanced Analysis Functions:}
#' \itemize{
#'   \item \code{\link{crossq.max}}: Compute maximum cross-quantilogram
#'   \item \code{\link{crossq.partial}}: Compute partial cross-quantilogram
#' }
#'
#' For a complete list of functions, see the package index.
#'
#' @references
#' Han, H., Linton, O., Oka, T., & Whang, Y. J. (2016). The cross-quantilogram:
#' Measuring quantile dependence and testing directional predictability between
#' time series. Journal of Econometrics, 193(1), 251-270.
#'
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom utils tail
#'
#' @docType package
#' @name quantilogram-package
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
