#' L2 Norm Tests for Matrix Effects in Non-Linear Calibration Curves
#'
#' This package provides statistical tests for detecting matrix effects in 
#' quadratic and 4-parameter logistic calibration curves using L2 norm distances.
#'
#' @docType package
#' @name MatrixEffects
NULL

# ---- Dependencies -----------------------------------------------------------
.require_namespace <- function(pkg){
  if(!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required.", call. = FALSE)
}

.require_namespace("sandwich")
.require_namespace("minpack.lm")

# ---- Small utilities --------------------------------------------------------

ws_params <- function(M, Sigma){
  e <- Re(eigen(M %*% Sigma, symmetric = FALSE, only.values = TRUE)$values)
  e <- e[e > 1e-12]; if(!length(e)) stop("All eigenvalues <= 0")
  tr1 <- sum(e); tr2 <- sum(e^2); list(lambda=e, df = tr1^2/tr2, c = tr1/(tr1^2/tr2))
}

apply_ridge <- function(A, ridge = 1e-8){
  if(any(!is.finite(A))) A[!is.finite(A)] <- 0
  if(rcond(A) < 1e-12) A <- A + diag(ridge * mean(diag(A)), nrow(A))
  A
}

# ---- Euclidean matrix effect calculation -----------------------------------

euclidean_matrix_effect <- function(fun_ref, fun_mat, x_range, center = FALSE) {
  tryCatch({
    L <- x_range[2] - x_range[1]
    
    if (center) {
      # Centered (baseline-invariant) version
      mean_diff <- integrate(function(x) fun_mat(x) - fun_ref(x), 
                             lower = x_range[1], upper = x_range[2])$value / L
      mean_ref <- integrate(function(x) fun_ref(x), 
                            lower = x_range[1], upper = x_range[2])$value / L
      
      diff_integral <- integrate(function(x) ((fun_mat(x) - fun_ref(x)) - mean_diff)^2,
                                 lower = x_range[1], upper = x_range[2])$value
      ref_integral <- integrate(function(x) (fun_ref(x) - mean_ref)^2,
                                lower = x_range[1], upper = x_range[2])$value
    } else {
      diff_integral <- integrate(function(x) (fun_mat(x) - fun_ref(x))^2,
                                 lower = x_range[1], upper = x_range[2])$value
      ref_integral <- integrate(function(x) fun_ref(x)^2,
                                lower = x_range[1], upper = x_range[2])$value
    }
    
    if (ref_integral < 1e-12) {
      warning("Reference curve has near-zero variation; delta2 undefined")
      return(NA_real_)
    }
    
    delta2 <- 100 * sqrt(diff_integral) / sqrt(ref_integral)
    return(delta2)
    
  }, error = function(e) {
    warning("Error in euclidean_matrix_effect: ", e$message)
    return(NA_real_)
  })
}

# ---- Quadratic model --------------------------------------------------------

gram_matrix_quad <- function(xmin, xmax, drop_intercept = TRUE, center = FALSE){
  M <- matrix(0,3,3); for(i in 0:2) for(j in 0:2){ p<-i+j+1; M[i+1,j+1]<-(xmax^p-xmin^p)/p }
  
  if(center){
    L <- xmax - xmin
    mu <- c(1, 
            (xmax^2 - xmin^2) / (2 * L),
            (xmax^3 - xmin^3) / (3 * L))
    M <- M - L * outer(mu, mu)
  }
  
  if(drop_intercept) M <- M[-1,-1]
  M
}

#' L2 Test for Quadratic Calibration Curves
#'
#' Performs an L2 norm test for matrix effects between quadratic calibration curves.
#'
#' @param x_ref Numeric vector of x values for reference curve
#' @param y_ref Numeric vector of y values for reference curve  
#' @param x_mat Numeric vector of x values for matrix curve
#' @param y_mat Numeric vector of y values for matrix curve
#' @param x_range Numeric vector of length 2 specifying integration range. If NULL, uses data range.
#' @param cov_type Character string specifying covariance estimator type
#' @param drop_intercept Logical indicating whether to drop intercept from comparison
#' @param center Logical indicating whether to use baseline-invariant (centered) test
#' @param alpha Significance level for the test
#'
#' @return An object of class "l2test" containing test results
#' @export
l2_test_quad <- function(x_ref, y_ref, x_mat, y_mat, 
                         x_range = NULL, 
                         cov_type = c("classical","const","HC0","HC1","HC2","HC3","HC4","HC4m","HC5"), 
                         drop_intercept = FALSE, 
                         center = FALSE, 
                         alpha = 0.05){
  
  cov_type <- match.arg(cov_type)
  
  # Fit quadratic models
  fit_ref <- lm(y_ref ~ x_ref + I(x_ref^2))
  fit_mat <- lm(y_mat ~ x_mat + I(x_mat^2))
  
  if(is.null(x_range)) x_range <- range(c(x_ref, x_mat))
  
  idx <- if(drop_intercept) 2:3 else 1:3
  delta <- coef(fit_ref)[idx] - coef(fit_mat)[idx]
  
  M <- gram_matrix_quad(x_range[1], x_range[2], drop_intercept, center)
  D <- as.numeric(t(delta) %*% M %*% delta)
  
  # Covariance estimation
  if(cov_type %in% c("classical", "const")) {
    Sig <- vcov(fit_ref)[idx,idx] + vcov(fit_mat)[idx,idx]
  } else {
    Sig <- sandwich::vcovHC(fit_ref, type = cov_type)[idx,idx] + 
      sandwich::vcovHC(fit_mat, type = cov_type)[idx,idx]
  }
  
  # Satterthwaite approximation + F test
  ws <- ws_params(M, Sig)
  df_den <- fit_ref$df.residual + fit_mat$df.residual
  F_stat <- D / (ws$c * ws$df)
  p <- 1 - pf(F_stat, ws$df, df_den)
  
  # Create prediction functions
  coef_ref <- coef(fit_ref)
  coef_mat <- coef(fit_mat)
  
  if (drop_intercept) {
    fun_ref <- function(x) coef_ref[2]*x + coef_ref[3]*x^2
    fun_mat <- function(x) coef_mat[2]*x + coef_mat[3]*x^2
  } else {
    fun_ref <- function(x) coef_ref[1] + coef_ref[2]*x + coef_ref[3]*x^2
    fun_mat <- function(x) coef_mat[1] + coef_mat[2]*x + coef_mat[3]*x^2
  }
  
  delta2 <- euclidean_matrix_effect(fun_ref, fun_mat, x_range, center = center)
  if (is.na(delta2)) delta2 <- 0
  
  structure(list(
    statistic = F_stat, p_value = p, D_statistic = D, df_num = ws$df, 
    df_denom = df_den, eigenvals = ws$lambda, scale_c = ws$c, 
    delta2 = delta2, x_range = x_range,
    fit_ref = fit_ref, fit_mat = fit_mat, fun_ref = fun_ref, fun_mat = fun_mat,
    method = sprintf("Quadratic L2 F-test (%s%s)", 
                     cov_type, 
                     if(center) ", centered" else ""), 
    alpha = alpha, model = "quad", center = center
  ), class = "l2test")
}

# ---- 4-PL functions ---------------------------------------------------------

f4 <- function(x,A,D,B,C) A + (D-A)/(1+(x/C)^(-B))

grad4 <- function(x,p){ A<-p[1];D<-p[2];B<-p[3];C<-p[4]; t<-(x/C)^(-B); d<-(1+t)^2; cbind(1-1/(1+t),1/(1+t),(D-A)*t*log(x/C)/d,-(D-A)*B*t/(C*d)) }

# Artificial-regression covariance -------------------------------------------

nlshc <- function(nlsfit, type = NULL){
  m <- nlsfit$m; X <- m$gradient(); r <- m$resid(); b <- coef(nlsfit)
  y_lin <- r + X %*% b; lm_fit <- lm(y_lin ~ X - 1)
  if(is.null(type) || type %in% c("classical", "const")) return(vcov(lm_fit))
  out <- sandwich::vcovHC(lm_fit, type = type)
  if(any(!is.finite(out))) stop(sprintf("vcovHC(%s) produced non-finite values", type))
  out
}


fit_4pl <- function(x, y, w=NULL, control=NULL){
  # Default control parameters
  default_control <- list(
    maxiter = 500,
    starts = list(
      c(min(y), max(y), 1, median(x)),
      c(quantile(y, 0.1), quantile(y, 0.9), 0.7, median(x)),
      c(min(y)*0.8, max(y)*1.2, 2, median(x))
    )
  )
  
  # Merge user control with defaults
  if(!is.null(control)) {
    for(name in names(control)) {
      default_control[[name]] <- control[[name]]
    }
  }
  control <- default_control
  
  for(st in control$starts){
    start_list <- as.list(st)
    names(start_list) <- c("A","D","B","C")
    
    fit <- try(minpack.lm::nlsLM(
      y ~ f4(x, A, D, B, C), 
      start = start_list,
      weights = w,
      control = minpack.lm::nls.lm.control(maxiter = control$maxiter)
    ), silent = TRUE)
    
    if(!inherits(fit, "try-error")) { 
      pm <- coef(fit); 
      if(all(pm[3:4] > 0)) return(fit) 
    }
  }
  stop("4-PL nlsLM failed")
}

#' L2 Test for 4-Parameter Logistic Calibration Curves
#'
#' Performs an L2 norm test for matrix effects between 4PL calibration curves.
#'
#' @param x_ref Numeric vector of x values for reference curve
#' @param y_ref Numeric vector of y values for reference curve  
#' @param x_mat Numeric vector of x values for matrix curve
#' @param y_mat Numeric vector of y values for matrix curve
#' @param x_range Numeric vector of length 2 specifying integration range.
#'   If NULL (default), uses overlap of data ranges. Use "union" for union of ranges.
#' @param n_grid Integer number of grid points for numerical integration (default 500)
#' @param measure Integration measure: "dose" for uniform in x, "logdose" for uniform in log(x)
#' @param cov_type Character string specifying covariance estimator type
#' @param center Logical indicating whether to use baseline-invariant (centered) test
#' @param ridge Numeric ridge parameter for matrix conditioning
#' @param alpha Significance level for the test
#' @param control List of control parameters for 4PL fitting
#'
#' @return An object of class "l2test" containing test results
#' @export
l2_test_4pl <- function(x_ref, y_ref, x_mat, y_mat, 
                        x_range = NULL,
                        n_grid = 500, 
                        measure = c("dose", "logdose"),
                        cov_type = c("classical","const","HC0","HC1","HC2","HC3","HC4","HC4m","HC5"), 
                        center = FALSE, 
                        ridge = 1e-8, 
                        alpha = 0.05, 
                        control = NULL){
  
  cov_type <- match.arg(cov_type)
  measure <- match.arg(measure)
  
  # Fit 4PL models
  fit_ref <- fit_4pl(x_ref, y_ref, rep(1, length(x_ref)), control)
  fit_mat <- fit_4pl(x_mat, y_mat, rep(1, length(x_mat)), control)
  
  pr <- coef(fit_ref); pm <- coef(fit_mat)
  
  # Determine integration range
  range_ref <- range(x_ref)
  range_mat <- range(x_mat)
  
  if (is.null(x_range)) {
    # Default: overlap (intersection) of ranges - avoids extrapolation
    x_range <- c(max(range_ref[1], range_mat[1]), 
                 min(range_ref[2], range_mat[2]))
    if (x_range[1] >= x_range[2]) {
      stop("Calibration ranges do not overlap. Specify x_range manually.")
    }
  } else if (identical(x_range, "union")) {
    x_range <- c(min(range_ref[1], range_mat[1]), 
                 max(range_ref[2], range_mat[2]))
  } else {
    if (length(x_range) != 2 || x_range[1] >= x_range[2]) {
      stop("x_range must be a numeric vector of length 2 with x_range[1] < x_range[2]")
    }
  }
  
  # Validate range for logdose measure
  if (measure == "logdose" && x_range[1] <= 0) {
    stop("x_range[1] must be positive for measure='logdose'")
  }
  
  # Generate integration grid
  if (measure == "dose") {
    x_grid <- seq(x_range[1], x_range[2], length.out = n_grid)
    dx <- (x_range[2] - x_range[1]) / (n_grid - 1)
    L <- x_range[2] - x_range[1]
  } else {
    x_grid <- exp(seq(log(x_range[1]), log(x_range[2]), length.out = n_grid))
    dx <- log(x_range[2] / x_range[1]) / (n_grid - 1)
    L <- log(x_range[2] / x_range[1])
  }
  
  # Compute curve difference
  d <- f4(x_grid, pr[1], pr[2], pr[3], pr[4]) - f4(x_grid, pm[1], pm[2], pm[3], pm[4])
  
  if (center) {
    mean_diff <- sum(d) * dx / L
    d <- d - mean_diff
  }
  
  Dobs <- sum(d*d) * dx
  
  # Compute Gram matrix
  J <- grad4(x_grid, (pr+pm)/2); M <- crossprod(J)*dx
  
  if (center) {
    mu_J <- colSums(J) * dx / L
    M <- M - L * tcrossprod(mu_J)
  }
  
  # Covariance estimation
  if(cov_type %in% c("classical", "const")) {
    V_ref <- nlshc(fit_ref, NULL)
    V_mat <- nlshc(fit_mat, NULL)
  } else {
    V_ref <- nlshc(fit_ref, cov_type)
    V_mat <- nlshc(fit_mat, cov_type)
  }
  
  Sigma <- apply_ridge(V_ref + V_mat, ridge)
  ws <- ws_params(M, Sigma)
  df_den <- length(x_ref) + length(x_mat) - 8
  F_stat <- Dobs/(ws$c*ws$df); p <- 1-pf(F_stat, ws$df, df_den)
  
  # Prediction functions
  fun_ref <- function(x) f4(x, pr[1], pr[2], pr[3], pr[4])
  fun_mat <- function(x) f4(x, pm[1], pm[2], pm[3], pm[4])
  
  # Calculate delta2 using same geometry
  y_ref_g <- fun_ref(x_grid)
  y_mat_g <- fun_mat(x_grid)
  
  if (!center) {
    num <- sqrt(sum((y_mat_g - y_ref_g)^2) * dx)
    den <- sqrt(sum(y_ref_g^2) * dx)
  } else {
    mean_diff_d2 <- sum(y_mat_g - y_ref_g) * dx / L
    mean_ref_d2 <- sum(y_ref_g) * dx / L
    
    diff_c <- (y_mat_g - y_ref_g) - mean_diff_d2
    ref_c <- y_ref_g - mean_ref_d2
    
    num <- sqrt(sum(diff_c^2) * dx)
    den <- sqrt(sum(ref_c^2) * dx)
  }
  
  delta2 <- if (den < 1e-12) NA_real_ else 100 * num / den
  if (is.na(delta2)) delta2 <- 0
  
  structure(list(
    statistic=F_stat, p_value=p, D_statistic=Dobs, df_num=ws$df, df_denom=df_den, 
    eigenvals=ws$lambda, scale_c=ws$c, delta2=delta2, x_range=x_range,
    fit_ref=fit_ref, fit_mat=fit_mat, fun_ref=fun_ref, fun_mat=fun_mat,
    method=sprintf("4-PL L2 F-test (%s%s, %s)", 
                   cov_type, 
                   if(center) ", centered" else "", 
                   measure), 
    alpha=alpha, model="4pl", center=center, measure=measure,
    ridge=ridge, n_grid=n_grid
  ), class="l2test")
}

# ---- Facade -----------------------------------------------------------------

#' L2 Test for Matrix Effects
#'
#' General interface for L2 norm tests that automatically detects model type.
#'
#' @param x_ref Numeric vector of x values for reference curve
#' @param y_ref Numeric vector of y values for reference curve  
#' @param x_mat Numeric vector of x values for matrix curve
#' @param y_mat Numeric vector of y values for matrix curve
#' @param model Character string specifying model type ("auto", "quad", "4pl")
#' @param ... Additional arguments passed to specific test functions
#'
#' @return An object of class "l2test" containing test results
#' @export
l2_test <- function(x_ref, y_ref, x_mat, y_mat, model=c("auto","quad","4pl"), ...){ 
  mdl <- match.arg(model)
  if(mdl == "auto") {
    x_range <- range(c(x_ref, x_mat))
    mdl <- if(x_range[2]/x_range[1] > 100) "4pl" else "quad"
  }
  switch(mdl, 
         quad = l2_test_quad(x_ref, y_ref, x_mat, y_mat, ...), 
         `4pl` = l2_test_4pl(x_ref, y_ref, x_mat, y_mat, ...)
  )
}

# ---- S3 Methods -------------------------------------------------------------

#' @export
print.l2test <- function(x, ...){ 
  cat(x$method, "\n")
  cat(sprintf("  F = %.4f  (df1 = %.2f, df2 = %.2f)  p = %.4g\n", x$statistic, x$df_num, x$df_denom, x$p_value))
  cat(sprintf("  D = %.4f\n", x$D_statistic))
  cat(sprintf("  delta2 (matrix effect) = %.2f%%\n", x$delta2))
  cat(sprintf("  X range: [%.3f, %.3f]\n", x$x_range[1], x$x_range[2]))
  invisible(x) 
}

#' Summary Method for L2 Test Results
#'
#' @param object An l2test object
#' @param ... Additional arguments
#'
#' @exportsummary.l2test <- function(object, ...) {
  cat("L2 Norm Test for Matrix Effects\n")
  cat("===============================\n\n")
  
  cat("Model:", toupper(object$model), "\n")
  cat("Method:", object$method, "\n")
  cat("X range:", sprintf("[%.3f, %.3f]", object$x_range[1], object$x_range[2]), "\n")
  if(!is.null(object$center) && object$center) {
    cat("Centering: Baseline-invariant (centered)\n")
  }
  cat("\n")
  
  cat("Test Results:\n")
  cat("  F-statistic:", sprintf("%.4f", object$statistic), "\n")
  cat("  Degrees of freedom:", sprintf("%.2f, %.2f", object$df_num, object$df_denom), "\n")
  sig <- symnum(object$p_value, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                symbols = c("***", "**", "*", ".", " "))
  cat("  p-value:", sprintf("%.4g", object$p_value), as.character(sig), "\n")
  cat("  Signif. codes: 0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  
  cat("Effect Size:\n")
  cat("  D-statistic:", sprintf("%.4f", object$D_statistic), "\n")
  cat("  delta2 (matrix effect):", sprintf("%.2f%%", object$delta2), "\n")
  cat("  Scale factor c:", sprintf("%.4g", object$scale_c), "\n\n")
  
  if(object$model == "quad") {
    cat("Quadratic Coefficients:\n")
    coef_ref <- coef(object$fit_ref)
    coef_mat <- coef(object$fit_mat)
    cat("  Reference: ", sprintf("%.6f + %.6f*x + %.6f*x^2", coef_ref[1], coef_ref[2], coef_ref[3]), "\n")
    cat("  Matrix:    ", sprintf("%.6f + %.6f*x + %.6f*x^2", coef_mat[1], coef_mat[2], coef_mat[3]), "\n")
  } else {
    cat("4PL Parameters:\n")
    pr <- coef(object$fit_ref); pm <- coef(object$fit_mat)
    cat("  Reference: A =", sprintf("%.4f", pr[1]), " D =", sprintf("%.4f", pr[2]), " B =", sprintf("%.4f", pr[3]), " C =", sprintf("%.4f", pr[4]), "\n")
    cat("  Matrix:    A =", sprintf("%.4f", pm[1]), " D =", sprintf("%.4f", pm[2]), " B =", sprintf("%.4f", pm[3]), " C =", sprintf("%.4f", pm[4]), "\n")
  }
  
  invisible(object)
}
