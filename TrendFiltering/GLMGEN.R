library(devtools)
library(testthat)
library(genlasso)
test_check("glmgen")
install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
set.seed(0)

trendfilter.control.list = function(rho=1, obj_tol=1e-5, obj_tol_newton=obj_tol,
                                    max_iter=200L, max_iter_newton=50L, 
                                    x_tol=1e-6, alpha_ls=0.5, gamma_ls=0.8,
                                    max_iter_ls=30L, tridiag=0) {
  
  z <- list(rho=rho, obj_tol=obj_tol, obj_tol_newton=obj_tol_newton,
            max_iter=max_iter, max_iter_newton=max_iter_newton, 
            x_tol=x_tol, alpha_ls=alpha_ls, gamma_ls=gamma_ls,
            max_iter_ls=max_iter_ls, tridiag=tridiag)
  z
}


trendfilter = function(x, y, weights, k = 2L,
                       family = c("gaussian", "logistic", "poisson"),
                       method = c("admm"),
                       beta0 = NULL,
                       lambda, nlambda = 50L, lambda.min.ratio = 1e-5,
                       thinning = NULL, verbose = F,
                       control = trendfilter.control.list(x_tol=1e-6*max(IQR(x),diff(range(x))/2))) {
  
  cl = match.call()
  family = match.arg(family)
  method = match.arg(method)
  family_cd = match(family, c("gaussian", "logistic", "poisson")) - 1L
  method_cd = match(method, c("admm")) - 1L
  
  if (missing(x) || is.null(x)) stop("x must be passed.")
  if (missing(y) || is.null(y)) { y = x; x = 1L:length(y) }
  else if (length(x) != length(y)) stop("x and y must have the same length.")
  n = length(y)
  ord = order(x)
  y = y[ord]
  x = x[ord]
  
  if (family_cd == 1L & !all(y %in% c(0,1)))
    warning("Logistic family should have all 0/1 responses.")
  if (family_cd == 2L & (any(round(y) != y) | any(y < 0)))
    warning("Poisson family requires non-negative integer responses.")
  
  if (missing(weights)) weights = rep(1L,length(y))
  if (any(weights==0)) stop("Cannot pass zero weights.")
  weights = weights[ord]
  
  if (is.na(family_cd)) stop("family argument must be one of 'gaussian', 'logistic', or 'poisson'.")
  if (k < 0 || k != floor(k)) stop("k must be a nonnegative integer.")
  if (n < k+2) stop("y must have length >= k+2 for kth order trend filtering.")
  if (k >= 3) warning("Large k leads to generally worse conditioning; k=0,1,2 are the most stable choices.")
  
  mindx = min(diff(x))
  if (!is.null(thinning) && !thinning && mindx == 0) {
    stop("Cannot pass duplicate x values; use observation weights, or use thinning=TRUE.")
  }
  
  # If the minimum difference between x points is < 1e-6 times
  # the interquartile range, then apply thinning, unless they
  # explicitly tell us not to
  if (mindx <= control$x_tol) {
    if (!is.null(thinning) && !thinning) {
      warning("The x values are ill-conditioned. Consider thinning. \nSee ?trendfilter for more info.")
    }
    else {	
      z = .Call("thin_R",
                sX = as.double(x),
                sY = as.double(y),
                sW = as.double(weights),
                sN = length(y),
                sK = as.integer(k),
                sControl = control,
                PACKAGE = "glmgen")
      x = z$x
      y = z$y
      weights = z$w
      n = z$n
      
      if (!is.null(beta0)) {
        z = .Call("thin_R",
                  sX = as.double(x),
                  sY = as.double(beta0),
                  sW = as.double(weights),
                  sN = length(y),
                  sK = as.integer(k),
                  sControl = control,
                  PACKAGE = "glmgen")
        beta0 = z$y
      }
    }
  }
  
  if (missing(lambda)) {
    if (nlambda < 1L || nlambda != floor(nlambda)) stop("nlambda must be a positive integer.")
    if (lambda.min.ratio <= 0 || lambda.min.ratio >= 1) stop("lamba.min.ratio must be between 0 and 1.")
    lambda = rep(0, nlambda)
    lambda_flag = FALSE
  } else {
    if (length(lambda) == 0L) stop("Must specify at least one lambda value.")
    if (min(lambda) < 0L) stop("All specified lambda values must be nonnegative.")
    if (any(order(lambda) != length(lambda):1L) & any(order(lambda) != 1L:length(lambda))) # ????
      warning("User-supplied lambda values should given in decending order for warm starts.")
    nlambda = length(lambda)
    lambda_flag = TRUE
  }
  if (!is.list(control) || (is.null(names(control)) && length(control) != 0L))
    stop("control must be a named list.")
  control = lapply(control, function(v) ifelse(is.numeric(v),
                                               as.double(v[[1]]), stop("Elements of control must be numeric.")))
  z = .Call("tf_R",
            sX = as.double(x),
            sY = as.double(y),
            sW = as.double(weights),
            sN = length(y),
            sK = as.integer(k),
            sFamily = as.integer(family_cd),
            sMethod = as.integer(method_cd),
            sBeta0 = beta0,
            sLamFlag = as.integer(lambda_flag),
            sLambda = as.double(lambda),
            sNlambda = as.integer(nlambda),
            sLambdaMinRatio = as.double(lambda.min.ratio),
            sVerbose = as.integer(verbose),
            sControl = control,
            PACKAGE = "glmgen")
  
  if (is.null(z)) stop("Unspecified error in C code.")
  colnames(z$beta) = as.character(round(z$lambda, 3))
  
  out = structure(list(y = y, x = x, weights = weights, k = as.integer(k),
                       lambda = z$lambda, beta0 = beta0, df = z$df, beta = z$beta, family = family,
                       method = method, n = length(y), p = length(y),
                       m = length(y) - as.integer(k) - 1L, obj = z$obj,
                       status = z$status, iter = z$iter, family=family, call = cl),
                  class = c("trendfilter","glmgen"))
  
  out
}

 n = 100
  x = runif(n, min=-2*pi, max=2*pi)
  y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
  

  xx = seq(min(x),max(x),length=100)
  lambda = out$lambda[25]
  yy = predict(out,x.new=xx,lambda=out$lambda[40])
  plot(x,y)
  lines(xx,yy,col=2)
  
  
  
  library(genlasso)
  library(glmgen)
  library(testthat)
  
  set.seed(0)
  EPS = 1e-4
  
  n = 100
  x = sort(runif(n, min=-2*pi, max=2*pi))
  y = 1.5*sin(x) + sin(2*x) + rnorm(n, sd=0.2)
  
  # Fused lasso w/o location values
  outGenlasso = genlasso::trendfilter(y, ord=0)
  outGlmgen = glmgen::trendfilter(y, k=0, lambda=outGenlasso$lambda)
  expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)
  
  # Fused lasso with location values
  outGenlasso = genlasso::trendfilter(y, x, ord=0)
  outGlmgen = glmgen::trendfilter(x, y, k=0, lambda=outGenlasso$lambda)
  expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)
  
  # Fused lasso using glmgen postions
  outGenlasso = genlasso::trendfilter(y, x, ord=0)
  outGlmgen = glmgen::trendfilter(x, y, k=0)
  predGenlasso = coef(outGenlasso, lambda=outGlmgen$lambda)$beta
  predGlmgen = predict(outGlmgen)
  expect_true(abs(max(predGenlasso - predGlmgen)) < EPS)
  
  p = trendfilter.control.list()
  p$max_iter = 2000
  p$obj_tol = 1e-10
  
  # Higher order trendfiltering w/o location values
  for (k in 1:2) {
    outGenlasso = genlasso::trendfilter(y, ord=k)
    outGlmgen = glmgen::trendfilter(y, k=k, lambda=outGenlasso$lambda, control=p)
    expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < EPS)
    print(abs(max(outGenlasso$beta - outGlmgen$beta)))
  }
  
  p$max_iter = 4000
  
  # Higher order trendfiltering with location values
  for (k in 1:2) {
    outGenlasso = genlasso::trendfilter(y, x, ord=k)
    outGlmgen = glmgen::trendfilter(x, y, k=k, lambda=outGenlasso$lambda, 
                                    control=p)
    expect_true(abs(max(outGenlasso$beta - outGlmgen$beta)) < 1e-2)
    print(abs(max(outGenlasso$beta - outGlmgen$beta)))
  }
  
  # Polynomial signal test
  lambda <- 1e4
  n = 4
  ks = c(0, 1, 2)
  for (k in ks) {
    X = seq.int(0, n)
    y = X^k
    mod <- trendfilter(x = X, y = y, k = k, lambda = lambda)
    out <- predict(mod, lambda = lambda, type = "response")
    expect_true(abs(max(out - y)) < EPS)
    print(abs(max(out - y)))
  }
  n=30
  X = 1:(n-1)
  beta = c(rep(0,n/2), rep (1,n/4), rep(0,n/4))
  y=1:(n-1)
  for (i in 1:(n)){
    y[i] = beta[i] + sqrt(X[i]/n)*rnorm(1,1)
  }
 
  