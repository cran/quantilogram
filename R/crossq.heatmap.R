##' Heatmap of Cross-Quantilogram
##'
##' This function creates a customizable heatmap visualization of the cross-quantilogram matrix
##' and returns a list containing the plot and a data frame of cross-quantilogram values with critical values.
##' The heatmap uses 0 values if the test of no correlation  
##' cannot be rejected, and it uses cross-quantilogram values otherwise.
##' The critical values are obtained by stationary bootstrap. 
##'
##' @param DATA An input matrix of dimensions T x 2, where T is the number of observations.
##'             Column 1 contains the first variable and Column 2 contains the second variable.
##'             This function will apply a k-period lag to the second variable during computation.
##' @param k An integer representing the lag.
##' @param vec.q A numeric vector of quantiles.
##' @param Bsize Bootstrap sample size for stationary bootstrap.
##' @param sigLev Significance level for statistical test. Default is 0.05 (5% significance level).
##' @param var1_name Name of the first variable (predicted variable). If NULL, defaults to "Variable 1".
##' @param var2_name Name of the second variable (predicting variable). If NULL, defaults to "Variable 2".
##' @param title Plot title. Default is "Cross-Quantilogram Heatmap".
##' @param subtitle Plot subtitle. Default is NULL (no subtitle).
##' @param colors A vector of colors for the heatmap. Default is c("blue", "lightblue", "white", "pink", "red").
##' @param color_values A vector of values for color scaling. Default is c(-1, -0.15, 0, 0.15, 1).
##' @param tile_border_color Color for tile borders. Default is "black".
##' @param tile_border_width Width for tile borders. Default is 0.5.
##' @param x_angle Angle for x-axis labels. Default is 90.
##' @param x_lab X-axis label. If NULL (default), it's automatically generated.
##' @param y_lab Y-axis label. If NULL (default), it's automatically generated.
##' @param legend_title Title for the legend. Default is "Cross-Q".
##'
##' @return A list containing two elements:
##'   \item{plot}{A ggplot object representing the cross-quantilogram heatmap.}
##'   \item{df.res}{A data frame containing cross-quantilogram values and critical values. It includes the following columns:
##'     \itemize{
##'       \item Quantile1: The quantile values for the first variable.
##'       \item Quantile2: The quantile values for the second variable.
##'       \item vCRQ: The cross-quantilogram values.
##'       \item Lower_CV: The lower critical values.
##'       \item Upper_CV: The upper critical values.
##'       \item Significant: A logical vector indicating whether the cross-quantilogram is significant at the given significance level.
##'     }
##'   }
##'
##' @import ggplot2
##' @importFrom scales rescale
##' @importFrom rlang .data
##' 
##' @references
##' Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016).
##' "The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series." \emph{Journal of Econometrics}, 193(1), 251-270.
##'
##' @examples
##' \dontrun{
##' ## data source 
##' data("sys.risk") 
##'
##' ## two variables data: T x 2 
##' DATA = sys.risk[,c("JPM", "Market")]
##' 
##' ## setup and estimation 
##' k = 1                             ## lag order 
##' vec.q  = seq(0.05, 0.95, 0.05)    ## a list of quantiles 
##' B.size = 200                      ## Repetition of bootstrap  
##' res = crossq.heatmap(DATA, k, vec.q, B.size) 
##' 
##' ## result 
##' print(res$plot)
##' }
##'
##' @author Heejoon Han, Oliver Linton, Tatsushi Oka and Yoon-Jae Whang
##' @export
##' 
crossq.heatmap = function(DATA, k, vec.q, Bsize, sigLev = 0.05,
                          var1_name = NULL, var2_name = NULL,
                          title = "Cross-Quantilogram Heatmap",
                          subtitle = NULL,
                          colors = c("blue", "lightblue", "white", "pink", "red"),
                          color_values = c(-1, -0.15, 0, 0.15, 1),
                          tile_border_color = "black",
                          tile_border_width = 0.5,
                          x_angle = 90,
                          x_lab = NULL,
                          y_lab = NULL,
                          legend_title = "Cross-Q") {
  
  if (missing(vec.q)) {
    stop("vec.q must be provided")
  }
  
  # Attempt to extract variable names from DATA if not provided
  if (is.null(var1_name) || is.null(var2_name)) {
    if (is.data.frame(DATA) && ncol(DATA) >= 2) {
      var1_name = ifelse(is.null(var1_name), names(DATA)[1], var1_name)
      var2_name = ifelse(is.null(var2_name), names(DATA)[2], var2_name)
    } else {
      var1_name = ifelse(is.null(var1_name), "Variable 1", var1_name)
      var2_name = ifelse(is.null(var2_name), "Variable 2", var2_name)
    }
  }
  
  # Calculate n.q once
  n.q = length(vec.q)
  
  # Pre-allocate memory for df.res
  df.res = data.frame(
    Quantile1   = numeric(n.q * n.q),
    Quantile2   = numeric(n.q * n.q),
    vCRQ        = numeric(n.q * n.q),
    Lower_CV    = numeric(n.q * n.q),
    Upper_CV    = numeric(n.q * n.q),
    Significant = logical(n.q * n.q)
  )

  
  # Counter for df.res
  counter = 1
  
  for (j1 in 1:n.q) {
    for (j2 in 1:n.q) {
      result = crossq.sb.opt(DATA, c(vec.q[j1], vec.q[j2]), k, Bsize, sigLev)
      
      # Determine significance
      is_significant = (result$vCRQ < result$vecCV[1]) | (result$vCRQ > result$vecCV[2])
      
      # Fill df.res
      df.res$Quantile1[counter]   = vec.q[j1]
      df.res$Quantile2[counter]   = vec.q[j2]
      df.res$vCRQ[counter]      = result$vCRQ
      df.res$Lower_CV[counter]    = result$vecCV[1]
      df.res$Upper_CV[counter]    = result$vecCV[2]
      df.res$Significant[counter] = is_significant
      
      counter = counter + 1
    }  
  }
  
  
  # Set default axis labels if not provided
  if (is.null(x_lab)) x_lab = paste("Quantile 2:", var2_name)
  if (is.null(y_lab)) y_lab = paste("Quantile 1:", var1_name)
  
  # Create the heatmap
  # - "y" for the 1st column variable and "x" for the 2nd column variable 
  p = ggplot(df.res, aes(y = .data$Quantile1, x = .data$Quantile2, 
                         fill = ifelse(.data$Significant, .data$vCRQ, 0))) +
    geom_tile(color = tile_border_color, linewidth = tile_border_width) +
    scale_fill_gradientn(
      colors = colors,
      values = scales::rescale(color_values),
      limits = c(-1, 1),
      breaks = seq(-1, 1, by = 0.2),
      name = legend_title
    ) +
    scale_x_continuous(breaks = 1:n.q,
                       labels = sprintf("%.2f", vec.q),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:n.q,
                       labels = sprintf("%.2f", vec.q),
                       expand = c(0, 0)) +
    theme_minimal() + 
    theme(
      axis.text.x = element_text(angle = x_angle, hjust = 1, vjust = 0.5),
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "right",
      axis.ticks = element_line(linewidth = 0.5),
      axis.ticks.length = unit(1.5, "mm"),
      axis.ticks.x.top = element_blank(),
      axis.ticks.y.right = element_blank(),
      legend.key.width = unit(1, "cm"),
      legend.text = element_text(hjust = 1)
    ) +
    coord_fixed() +
    labs(y = y_lab, x = x_lab, title = title, subtitle = subtitle)
  
  # Return both the plot and the data
  return(list(plot = p, df.res = df.res))
}




