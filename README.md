# MatrixEffects

L2 norm tests for matrix effects in non-linear calibration curves.

## Overview

MatrixEffects provides statistical tests for detecting matrix effects in quadratic and four-parameter logistic (4PL) calibration curves using L2 norm distances.

## Installation

You can install the development version from GitHub:

```r
# Install devtools
install.packages("devtools")

# Install MatrixEffects
devtools::install_github("BrevikM/MatrixEffects")
```

## Quick Start

```r
library(MatrixEffects)

# Example data
x_ref <- c(1, 20, 40, 60, 80, 100)
y_ref <- c(0.1, 0.8, 2.1, 4.2, 7.1, 10.8)
x_mat <- c(1, 20, 40, 60, 80, 100)  
y_mat <- c(0.1, 0.9, 2.4, 4.8, 8.1, 12.2)

# Perform the L2 F-test for quadratic calibration
result <- l2_test(x_ref, y_ref, x_mat, y_mat, model = "quad")
print(result)

# View detailed results
summary(result)
```

## Key Features

### Robust Statistical Testing
- F-distributed test statistics for both quadratic and 4PL models
- Multiple heteroscedasticity-consistent covariance estimators (HC0-HC5)
- Proper handling of non-linear parameter uncertainties

### Unified Effect Size Metric
- δ₂ (delta-2): Generalized matrix effect metric
- Reduces to classical slope ratio for linear calibrations
- Captures overall shape differences for non-linear curves

###
