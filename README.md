# Cross-Quantilogram

## Description

The `quantilogram` package provides estimation and inference methods for the cross-quantilogram. The cross-quantilogram is a measure of nonlinear dependence between two variables, based on either unconditional or conditional quantile functions. It can be considered an extension of the correlogram, which is a correlation function over multiple lag periods that mainly focuses on linear dependency.

This package allows users to detect the presence of directional predictability from one time series to another and provides a statistical inference method based on the stationary bootstrap.

## Installation

You can install the released version of quantilogram from [CRAN](https://CRAN.R-project.org) with:
  
  ```R
install.packages("quantilogram")
```

## Usage

Here's a basic example of how to use the quantilogram package:

```R
library(quantilogram)

# Load example data
data("sys.risk")

# Select two variables
D = sys.risk[, c("JPM", "Market")]

# Set parameters
k = 1                             # lag order 
vec.q = seq(0.05, 0.95, 0.05)     # a list of quantiles 
B.size = 200                      # Repetition of bootstrap  

# Compute and plot cross-quantilogram
res = heatmap.crossq(D, k, vec.q, B.size) 

# Display the plot
print(res$plot)
```

For more detailed examples and function descriptions, please refer to the package documentation.

## References

The methods implemented in this package are based on the following key publications:

1. Linton, O., and Whang, Y. J. (2007). The quantilogram: With an application to evaluating directional predictability. Journal of Econometrics, 141(1), 250-282. [doi:10.1016/j.jeconom.2007.01.004](https://doi.org/10.1016/j.jeconom.2007.01.004)

2. Han, H., Linton, O., Oka, T., and Whang, Y. J. (2016). The cross-quantilogram: Measuring quantile dependence and testing directional predictability between time series. Journal of Econometrics, 193(1), 251-270. [doi:10.1016/j.jeconom.2016.03.001](https://doi.org/10.1016/j.jeconom.2016.03.001)

## License

This package is free and open source software, licensed under GPL (>= 3).
