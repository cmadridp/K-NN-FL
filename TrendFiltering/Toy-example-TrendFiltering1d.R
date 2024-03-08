rm(list=ls())
library(devtools)
library(testthat)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(extrafont)
library(tidyverse)
library(GGally)
library(fields)
library(dplyr)
library(tidyr)
library(faux)
library(mvtnorm)
library(datasets)
library(gridExtra)  
library(R.matlab)
library(vegan) 
#GenLasso
library(devtools)
library(testthat)
library(genlasso)

test_check("glmgen")
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

###Cross validation
cv=function(X_s_training,y_i_js_training,folds_number,new_lambda){
  ms_train_test<-matrix(NaN,nrow=n, ncol=1)
  for(i in 1:(n)){
    ms_train_test[i]=(ms_training[i])/(folds_number)
  }
  nm_train_test=sum(ms_train_test)
  
  ms_train_train<-matrix(NaN,nrow=n, ncol=1)
  for(i in 1:(n)){
    ms_train_train[i]=4*(ms_training[i])/(folds_number)
  }
  nm_train_train=sum(ms_train_train)
  mse_lambdas=matrix(NaN,folds_number,length(new_lambda))
  lambda_min=matrix(NaN,1,length(new_lambda))
  for(k in 1:folds_number){
    X_s_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    y_i_js_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    Weights_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    
    X_s_training_test=matrix(NaN,nrow=nm_train_test,ncol=1)
    y_i_js_training_test=matrix(NaN,nrow=nm_train_test,ncol=1)
    Weights_training_test=matrix(NaN,nrow=nm_train_test,ncol=1)
    aux_2=0
    aux_1=0
    aux_3=0
    aux_4=0
    aux_5=0
    aux_6=0
    for (i in 1:n) {
      #print(i)
      aux=(ms_training[i])/(folds_number)*(k-1)
      aux_0=(ms_training[i])/(folds_number)*(k-1)
      if(aux>0){
        X_s_training_training[(aux_2+1):(aux_2+aux_0),]=X_s_training[(aux_1+1):(aux_1+(k-1)*(ms_training[i])/(folds_number)),]
        y_i_js_training_training[(aux_2+1):(aux_2+aux_0),]=y_i_js_training[(aux_1+1):(aux_1+(k-1)*(ms_training[i])/(folds_number)),]
        Weights_training_training[(aux_2+1):(aux_2+aux_0),]=Weights_training[(aux_1+1):(aux_1+(k-1)*(ms_training[i])/(folds_number)),]
        
        X_s_training_test[(aux_5+1):((ms_training[i])/(folds_number)+aux_5),]=X_s_training[(aux_4+(k-1)*(ms_training[i])/(folds_number)+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        y_i_js_training_test[(aux_5+1):((ms_training[i])/(folds_number)+aux_5),]=y_i_js_training[(aux_4+(k-1)*(ms_training[i])/(folds_number)+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        Weights_training_test[(aux_5+1):((ms_training[i])/(folds_number)+aux_5),]=Weights_training[(aux_4+(k-1)*(ms_training[i])/(folds_number)+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        
        if(aux+(ms_training[i])/(folds_number)+1<=ms_training[i]){
          X_s_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=X_s_training[(aux_4+k*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
          y_i_js_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=y_i_js_training[(aux_4+k*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
          Weights_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=Weights_training[(aux_4+k*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        }
        aux_2=aux_2+ms_train_train[i]
        aux_1=aux_1+ms_training[i]
        aux_5=aux_5+ms_train_test[i]
      }
      else{
        X_s_training_test[(aux_5+aux_0+1):((k)*(ms_training[i])/(folds_number)+aux_5),]=X_s_training[(aux_4+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        y_i_js_training_test[(aux_5+aux_0+1):((k)*(ms_training[i])/(folds_number)+aux_5),]=y_i_js_training[(aux_4+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        Weights_training_test[(aux_5+aux_0+1):((k)*(ms_training[i])/(folds_number)+aux_5),]=Weights_training[(aux_4+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        
        
        X_s_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3)]=X_s_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        y_i_js_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3)]=y_i_js_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        Weights_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3)]=Weights_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        aux_2=aux_2+ms_train_train[i]
        aux_5=aux_5+ms_train_test[i]
        aux_1=aux_1+ms_training[i]
      }
      #aux=ms_training[i]+aux
      aux_4=aux_4+ms_training[i]
      aux_3=aux_3+ms_train_train[i]
      aux_6=aux_6+ms_train_test[i]
    }
    #creating lambdas
    #lambda=seq(0.00001,2,0.001)
    for(i in 1:(length(new_lambda))){
      out = trendfilter(X_s_training_training, y_i_js_training_training, Weights_training_training, k = 0,lambda = new_lambda)
      xx = X_s_training_test
      yy = predict(out,x.new=xx,lambda=out$lambda[i])
      mse_lambdas[k,i]=MSE(yy,y_i_js_training_test)
    }
    
    #plot(X_s_training_test,y_i_js_training_test,xlim=c(0, 1),ylim=c(0,10))
    #lines(xx,yy,col=2,lwd=3)
  }
  for(i in 1:(length(new_lambda))){
    lambda_min[i]=sum(mse_lambdas[,i])/folds_number
  }
  lambda_opt=which(lambda_min == min(lambda_min))
  return(lambda_opt)
}

###Cross validation
cv_smooth_spline=function(X_s_training,y_i_js_training,folds_number,new_lambda){
  ms_train_test<-matrix(NaN,nrow=n, ncol=1)
  for(i in 1:(n)){
    ms_train_test[i]=(ms_training[i])/(folds_number)
  }
  nm_train_test=sum(ms_train_test)
  
  ms_train_train<-matrix(NaN,nrow=n, ncol=1)
  for(i in 1:(n)){
    ms_train_train[i]=4*(ms_training[i])/(folds_number)
  }
  nm_train_train=sum(ms_train_train)
  mse_lambdas=matrix(NaN,folds_number,length(new_lambda))
  lambda_min=matrix(NaN,1,length(new_lambda))
  for(k in 1:folds_number){
    X_s_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    y_i_js_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    Weights_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    
    X_s_training_test=matrix(NaN,nrow=nm_train_test,ncol=1)
    y_i_js_training_test=matrix(NaN,nrow=nm_train_test,ncol=1)
    Weights_training_test=matrix(NaN,nrow=nm_train_test,ncol=1)
    aux_2=0
    aux_1=0
    aux_3=0
    aux_4=0
    aux_5=0
    aux_6=0
    for (i in 1:n) {
      #print(i)
      aux=(ms_training[i])/(folds_number)*(k-1)
      aux_0=(ms_training[i])/(folds_number)*(k-1)
      if(aux>0){
        X_s_training_training[(aux_2+1):(aux_2+aux_0),]=X_s_training[(aux_1+1):(aux_1+(k-1)*(ms_training[i])/(folds_number)),]
        y_i_js_training_training[(aux_2+1):(aux_2+aux_0),]=y_i_js_training[(aux_1+1):(aux_1+(k-1)*(ms_training[i])/(folds_number)),]
        Weights_training_training[(aux_2+1):(aux_2+aux_0),]=Weights_training[(aux_1+1):(aux_1+(k-1)*(ms_training[i])/(folds_number)),]
        
        X_s_training_test[(aux_5+1):((ms_training[i])/(folds_number)+aux_5),]=X_s_training[(aux_4+(k-1)*(ms_training[i])/(folds_number)+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        y_i_js_training_test[(aux_5+1):((ms_training[i])/(folds_number)+aux_5),]=y_i_js_training[(aux_4+(k-1)*(ms_training[i])/(folds_number)+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        Weights_training_test[(aux_5+1):((ms_training[i])/(folds_number)+aux_5),]=Weights_training[(aux_4+(k-1)*(ms_training[i])/(folds_number)+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        
        if(aux+(ms_training[i])/(folds_number)+1<=ms_training[i]){
          X_s_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=X_s_training[(aux_4+k*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
          y_i_js_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=y_i_js_training[(aux_4+k*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
          Weights_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=Weights_training[(aux_4+k*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        }
        aux_2=aux_2+ms_train_train[i]
        aux_1=aux_1+ms_training[i]
        aux_5=aux_5+ms_train_test[i]
      }
      else{
        X_s_training_test[(aux_5+aux_0+1):((k)*(ms_training[i])/(folds_number)+aux_5),]=X_s_training[(aux_4+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        y_i_js_training_test[(aux_5+aux_0+1):((k)*(ms_training[i])/(folds_number)+aux_5),]=y_i_js_training[(aux_4+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        Weights_training_test[(aux_5+aux_0+1):((k)*(ms_training[i])/(folds_number)+aux_5),]=Weights_training[(aux_4+1):(aux_4+(k)*(ms_training[i])/(folds_number)),]
        
        
        X_s_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3)]=X_s_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        y_i_js_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3)]=y_i_js_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        Weights_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3)]=Weights_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
        aux_2=aux_2+ms_train_train[i]
        aux_5=aux_5+ms_train_test[i]
        aux_1=aux_1+ms_training[i]
      }
      #aux=ms_training[i]+aux
      aux_4=aux_4+ms_training[i]
      aux_3=aux_3+ms_train_train[i]
      aux_6=aux_6+ms_train_test[i]
    }
    #creating lambdas
    #lambda=seq(0.00001,2,0.001)
    for(i in 1:(length(new_lambda))){
      out = smooth.spline(X_s_training_training, y_i_js_training_training, Weights_training_training,lambda = new_lambda[i])
      xx = X_s_training_test
      yy = predict(out,xx)
      yy_1=yy$y
      mse_lambdas[k,i]=MSE(yy_1,y_i_js_training_test)
    }
    
    #plot(X_s_training_test,y_i_js_training_test,xlim=c(0, 1),ylim=c(0,10))
    #lines(xx,yy,col=2,lwd=3)
  }
  for(i in 1:(length(new_lambda))){
    lambda_min[i]=sum(mse_lambdas[,i])/folds_number
  }
  lambda_opt=which(lambda_min == min(lambda_min))
  return(lambda_opt)
}



MSE=function(y_test,yy){
  res=0
  for (i in 1:(length(yy))) {
    res=res+(yy[i]-y_test[i])^2
  }
  return(res/length(yy))
}
mse_new <- function(y_predicted, X_s_test, Weights_tes,constantb) {
  diff <- y_predicted[,1] - f_vec_eva(X_s_test,constantb)
  mse <- sum((diff^2) * Weights_tes[,1])
  return(mse)
}
#############################
#############################
######### Scenario 1 ######## Piecewise constant
#############################
#############################
f_vec_eva <- function(x,constantb) {
  d <- ncol(x)
  n <- nrow(x)
  result <- rep(0, n)
  ones_1 <- rep(1, d)
  for (i in 1:n) {
    result[i] <- f(x[i, ],constantb)
  }
  return(result)
}

f_aux=function(x){
  if(x<=0.2){
    return(5)
  }
  else if(0.2<x && x<=0.4){
    return(2)
  }
  else if(0.4<x && x<=0.6){
    return(8)
  }
  else{
    return(1)
  }
}

#The function choice
f=function(x,b){
  if(x<=0.2){
    return(5-b)
  }
  else if(0.2<x && x<=0.4){
    return(2-b)
  }
  else if(0.4<x && x<=0.6){
    return(8-b)
  }
  else{
    return(1-b)
  }
}

#we set the parameters
r=0
n=200
d=1
mult=2
m1=10*(mult+2)
m2=10*(mult)
m3=10*(mult+3)
m4=10*(mult+1)
#we define the n_ts
ms<-matrix(NaN,nrow=n, ncol=1)
for(i in 1:(n/4+1)){
  ms[i]=m1
}
for(i in (n/4+1):(n/2+1)){
  ms[i]=m2
}
for(i in (n/2+1):(3*n/4+1)){
  ms[i]=m3
}
for(i in (3*n/4+1):n){
  ms[i]=m4
}
nm=sum(ms)



r=0
repetitions=1
result_rep=matrix(NaN,repetitions,16)
result_rep_smoothspline=matrix(NaN,repetitions,16)
start_time <- Sys.time()
mult1=c(0.5,1,1.5,2)
n1=c(200,400,600,800)
for (s in 1:4) {
  
  n=n1[s]
  for (l in 1:4) {
    r=r+1
    #print(r)
    #print(l)
    mult=mult1[l]
    m1=10*(mult+2)
    m2=10*(mult)
    m3=10*(mult+3)
    m4=10*(mult+1)
    #we define the n_ts
    ms<-matrix(NaN,nrow=n, ncol=1)
    for(i in 1:(n/4+1)){
      ms[i]=m1
    }
    for(i in (n/4+1):(n/2+1)){
      ms[i]=m2
    }
    for(i in (n/2+1):(3*n/4+1)){
      ms[i]=m3
    }
    for(i in (3*n/4+1):n){
      ms[i]=m4
    }
    nm=sum(ms)
    for(rep in 1:repetitions){
      set.seed(rep)
      X_s<-matrix(NaN,nrow=nm, ncol=d)
      set.seed(rep)
      for(i in 1:nm){
        for(j in 1:d)
        {
          X_s[i,j]=runif(1,0,1)
        }
      }
      X_s
      ##########empirical mean=0
      constantb=0
      for (i in 1:nm) {
        for (j in 1:d) {
          constantb=constantb+f_aux(X_s[i,j])/nm
        }
      }
      
      ##creating the Delta_i functions
      h_t<-function(x_i_j,t){
        aux=1
        for(i in 1:(d))
        {
          aux=sqrt(2)*sin(t*(pi/2)*x_i_j[i])*aux
        }
        return(aux)
      }
      #running h_t
      aux=h_t(0.2,5)
      
      #defining the b_t_is
      b_i_ts=matrix(NaN,nrow=25,ncol=nm)
      #b_t_is=rnorm(btjs_val,0,1)
      for(i in 1:nm)
      {
        for (t in 1:25) {
          b_i_ts[t,i]=rnorm(1,0,1)
        }
      }
      b_i_ts
      #defining the \xi_t
      xi_i=function(x,i){
        aux=0
        #defining the b_t_is
        #b_t_js=rnorm(50,0,1)
        for (t in 1:25) {
          aux=aux+(b_i_ts[t,i]/(t^2))*h_t(x,t)
        }
        return(aux)
      }
      #we create the epsilon_ij--------------------------------------
      epsilon_s<-matrix(NaN,nrow=nm, ncol=1)
      epsilon_s=rnorm(nm,0,1)
      
      #we create the y_i_j----------------------------------------
      y_i_js=matrix(NaN,nrow=nm,ncol=1)
      Weights=matrix(NaN,nrow=nm,ncol=1)
      aux=0
      for(i in 1:n)
      {
        aux=ms[i]+aux
        for(j in 1:(ms[i]))
        {
          Weights[aux-ms[i]+j]=1/(n*ms[i])
          y_i_js[aux-ms[i]+j]=f(X_s[aux-ms[i]+j,],constantb)+xi_i(X_s[aux-ms[i]+j,],i)+epsilon_s[aux-ms[i]+j]
        }
      }
      y_i_js
      Weights
      
      
      
      
      #######################################
      #######################################
      ######creating data_training ##########
      #######################################
      #######################################
      ms_training<-matrix(NaN,nrow=n, ncol=1)
      for(i in 1:n){
        ms_training[i]=ms[i]-5*(ms[i]/10-mult+1)
      }
      nm_aux=sum(ms_training)
      
      
      ##### create X_s for training
      X_s_training=matrix(NaN,nrow=nm_aux,ncol=1)
      aux=0
      aux_1=0
      for(i in 1:n)
      {
        X_s_training[(aux+1):(aux+ms_training[i]),]=X_s[(aux_1+1):(aux_1+ms_training[i]),]
        aux_1=ms[i]+aux_1
        aux=ms_training[i]+aux
      }
      #we create the y_i_j for training----------------------------------------
      y_i_js_training=matrix(NaN,nrow=nm_aux,ncol=1)
      Weights_training=matrix(NaN,nrow=nm_aux,ncol=1)
      aux=0
      aux_1=0
      for(i in 1:n)
      {
        Weights_training[(aux+1):(aux+ms_training[i]),]=Weights[(aux_1+1):(aux_1+ms_training[i]),]
        y_i_js_training[(aux+1):(aux+ms_training[i]),]=y_i_js[(aux_1+1):(aux_1+ms_training[i]),]
        aux_1=ms[i]+aux_1
        aux=ms_training[i]+aux
      }
      
      #######################################
      #######################################
      ######creating data_test ##############
      #######################################
      #######################################
      ms_test<-matrix(NaN,nrow=n, ncol=1)
      for(i in 1:n){
        ms_test[i]=5*(ms[i]/10-mult+1)
      }
      nm_aux_test=sum(ms_test)
      
      
      ##### create X_s for test
      X_s_test=matrix(NaN,nrow=nm_aux_test,ncol=1)
      aux=0
      aux_1=0
      for(i in 1:n)
      {
        aux_1=ms[i]-ms_test[i]+aux_1
        X_s_test[(aux+1):(aux+ms_test[i]),]=X_s[(aux_1+1):(aux_1+ms_test[i]),]
        aux_1=aux_1+ms_test[i]
        aux=ms_test[i]+aux
      }
      #we create the y_i_j for test----------------------------------------
      y_i_js_test=matrix(NaN,nrow=nm_aux_test,ncol=1)
      Weights_test=matrix(NaN,nrow=nm_aux_test,ncol=1)
      aux=0
      aux_1=0
      for(i in 1:n)
      {
        aux_1=ms[i]-ms_test[i]+aux_1
        Weights_test[(aux+1):(aux+ms_test[i]),]=Weights[(aux_1+1):(aux_1+ms_test[i]),]
        y_i_js_test[(aux+1):(aux+ms_test[i]),]=y_i_js[(aux_1+1):(aux_1+ms_test[i]),]
        aux_1=aux_1+ms_test[i]
        aux=ms_test[i]+aux
      }
      
      #######################################
      #######################################
      ###### Executing Filtering ############
      #######################################
      #######################################
      folds_number=5
      lambda=c( 5.506572*exp(-7),exp(-7)
      )
      new_lambda= sort(lambda, decreasing=TRUE)
      lambda_opt=cv(X_s_training,y_i_js_training,folds_number,new_lambda)
      out = trendfilter(X_s_training, y_i_js_training, Weights_training, k = 0,lambda = new_lambda[lambda_opt])
      xx = X_s_test
      yy = predict(out,x.new=xx,lambda=new_lambda[lambda_opt])
      result_rep[rep,r]=mse_new(yy,X_s_test,Weights_test,constantb)
      
      
      lambda_smooth_spline=seq(0.00000000001,0.00000005,0.000005)
      lambda_opt_spline=cv_smooth_spline(X_s_training,y_i_js_training,folds_number,lambda_smooth_spline)
      ou_spline=smooth.spline(X_s_training, y_i_js_training, w=Weights_training,lambda=lambda_smooth_spline[lambda_opt_spline])
      xx = X_s_test
      yy = predict(ou_spline,xx)
      yy_1=yy$y
      result_rep_smoothspline[rep,r]=mse_new(yy_1,X_s_test,Weights_test,constantb)
      print(rep)
    }}}
end_time <- Sys.time()
end_time - start_time

result_rep
result_rep_smoothspline
