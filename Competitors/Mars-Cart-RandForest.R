

rm(list=ls())
library(randomForest)
library(vegan)
#library(rsample)   # data splitting
library(earth)     # fit MARS models
#library(caret)     # automating the tuning process
library(vip)       # variable importance
library(pdp)
library(mda)
library(rpart)
# variable relationships
###Cross validation
cv_mars=function(X_s_training,y_i_js_training,folds_number,number_of_neighb,penalty){
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
  new_penalty=sort(penalty)
  mse_penalty=matrix(NaN,folds_number,length(new_penalty))
  penalty_min=matrix(NaN,1,length(new_penalty))
  for(k in 1:folds_number){
    X_s_training_training=matrix(NaN,nrow=nm_train_train,ncol=d)
    y_i_js_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    Weights_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    
    X_s_training_test=matrix(NaN,nrow=nm_train_test,ncol=d)
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
        
        
        X_s_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=X_s_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
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
    for(i in 1:(length(new_penalty))){
      fit=mars(X_s_training_training, y_i_js_training_training, Weights=Weights_training_training,penalty=new_penalty[i])
      yy=predict(fit,X_s_training_test,penalty=new_penalty[i])
      mse_penalty[k,i]=mse_new(yy,X_s_training_test,Weights_training_test)
    }
  }
  for(i in 1:(length(new_penalty))){
    penalty_min[i]=sum(mse_penalty[,i])/folds_number
  }
  penalty_opt=which(penalty_min == min(penalty_min))
  return(penalty_opt)
}

cv_rpart=function(X_s_training,y_i_js_training,folds_number,number_of_neighb,cp){
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
  new_cp=sort(cp)
  mse_cp=matrix(NaN,folds_number,length(new_cp))
  cp_min=matrix(NaN,1,length(new_cp))
  for(k in 1:folds_number){
    X_s_training_training=matrix(NaN,nrow=nm_train_train,ncol=d)
    y_i_js_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    Weights_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    
    X_s_training_test=matrix(NaN,nrow=nm_train_test,ncol=d)
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
        
        
        X_s_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=X_s_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
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
    for(i in 1:(length(new_cp))){
      fit1= rpart(y_i_js_training_training~., data = data.frame(X_s_training_training, y_i_js_training_training), weights=Weights_training_training,
                  control = rpart.control(new_cp[i]))
      fit.pruned = prune(fit1, cp = new_cp[i])
      pred_y = predict(fit.pruned, data.frame(X_s_training_test))
      pred_y=as.matrix(pred_y)
      mse_cp[k,i]=mse_new(pred_y,X_s_training_test,Weights_training_test)
    }
  }
  for(i in 1:(length(new_cp))){
    cp_min[i]=sum(mse_cp[,i])/folds_number
  }
  cp_opt=which(cp_min == min(cp_min))
  return(cp_opt)
}


cv_randomForest=function(X_s_training,y_i_js_training,folds_number,number_of_neighb,nodesize){
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
  new_nodesize=sort(nodesize)
  mse_nodesize=matrix(NaN,folds_number,length(new_nodesize))
  nodesize_min=matrix(NaN,1,length(new_nodesize))
  for(k in 1:folds_number){
    X_s_training_training=matrix(NaN,nrow=nm_train_train,ncol=d)
    y_i_js_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    Weights_training_training=matrix(NaN,nrow=nm_train_train,ncol=1)
    
    X_s_training_test=matrix(NaN,nrow=nm_train_test,ncol=d)
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
        
        
        X_s_training_training[(aux_2+aux_0+1):(ms_train_train[i]+aux_3),]=X_s_training[(aux_4+(k)*(ms_training[i])/(folds_number)+1):(ms_training[i]+aux_4),]
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
    for(i in 1:(length(new_nodesize))){
      fit2=randomForest(x=X_s_training_training,y=y_i_js_training_training,Weights=Weights_training_training,ntree = 500,nodesize=nodesize[i])
      predy1=predict(fit2,newdata =X_s_training_test,type="response",ntree = 500,nodesize=nodesize[i])
      predy1=as.matrix(predy1)
      mse_nodesize[k,i]=mse_new(predy1,X_s_training_test,Weights_training_test)
    }
  }
  for(i in 1:(length(new_nodesize))){
    nodesize_min[i]=sum(mse_nodesize[,i])/folds_number
  }
  nodesize_opt=which(nodesize_min == min(nodesize_min))
  return(nodesize_opt)
}






MSE=function(y_test,yy){
  res=0
  for (i in 1:(length(yy))) {
    res=res+(yy[i]-y_test[i])^2
  }
  return(res/length(yy))
}

mse_new <- function(y_predicted, X_s_test, Weights_tes) {
  diff <- y_predicted[,1] - f_vec_eva(X_s_test)
  mse <- sum((diff^2) * Weights_tes)
  return(mse)
}






nearest_neighbors = function(x,obs, k, FUN, p = NULL){
  
  # Check the number of observations is the same
  if(ncol(x) != ncol(obs)){
    stop('Data must have the same number of variables')
  }
  
  # Calculate distance, considering p for Minkowski
  if(is.null(p)){
    dist = apply(x,1, FUN,obs)
  }else{
    dist = apply(x,1, FUN,obs,p)
  }
  
  # Find closest neighbours
  distances = sort(dist)[1:k]
  neighbor_ind = which(dist %in% sort(dist)[1:k])
  
  if(length(neighbor_ind)!= k){
    warning(
      paste('Several variables with equal distance. Used k:',length(neighbor_ind))
    )
  }
  
  ret = list(neighbor_ind, distances)
  return(ret)
}
euclidean_distance = function(a, b){
  #  We check that they have the same number of observation
  if(length(a) == length(b)){
    sqrt(sum((a-b)^2))
  } else{
    stop('Vectors must be of the same length')
  }
}


K_function=function(x_i,x,X_ij_training){
  aux=as.matrix(x)
  aux=t(aux)
  ind=nearest_neighbors(X_ij_training,aux,5,euclidean_distance)[[1]]
  for (j in 1:5) {
    if(X_ij_training[ind[j],1]==x_i[1] && X_ij_training[ind[j],2]==x_i[2]){
      return(1)
    }
  }
  return(0)
}

predict_function_for_x=function(x,X_ij_training,nm_training,fitted_with_training){
  aux=0
  for (i in 1:nm_training) {
    aux=aux+K_function(X_ij_training[i,],x,X_ij_training)
  }
  aux_1=0
  for (i in 1:nm_training) {
    aux_1=aux_1+fitted_with_training$fit[i]*K_function(X_ij_training[i,],x,X_ij_training)/aux
  }
  return(aux_1)
}
predict_function_for_a_set=function(X_test,nm_test,X_ij_training,nm_training,fitted_with_training){
  ans=matrix(NaN,nm_test,1)
  for(i in 1:nm_test){
    ans[i]=predict_function_for_x(X_test[i,],X_ij_training,nm_training,fitted_with_training)
  }
  return(ans)
}



#############################
#############################
######### Scenario 1 ######## Piecewise constant
#############################
#############################




######The constants
f <- function(x) {
  # Ensure x is treated as a matrix
  if (is.vector(x)) {
    x <- matrix(x, nrow = 1)
  }
  
  # Determine the number of columns in x
  d <- ncol(x)
  
  # Create a row vector of ones with the same number of columns as x
  ones_1 <- rep(1, d)
  
  # Calculate the squared Euclidean norms
  norm1 <- sum((x - 1/3 * ones_1)^2)
  norm2 <- sum((2 * x - 5/3 * ones_1)^2)
  
  # Compare the two norms and return the result
  if (norm1 < norm2) {
    return(1)
  } else {
    return(-1)
  }
}



f_vec_eva <- function(x) {
  d <- ncol(x)
  n <- nrow(x)
  result <- rep(0, n)
  ones_1 <- rep(1, d)
  for (i in 1:n) {
    result[i] <- f(x[i, ])
  }
  return(result)
}




#we set the parameters
r=0
n=1000
d=2
mult=0.5
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



####################
#################### SIMULATION WITH REPETITION
####################
r=0
d=2
repetitions=1
result_rep_mars=matrix(NaN,repetitions,1)
result_rep_cart=matrix(NaN,repetitions,1)
result_rep_RandFor=matrix(NaN,repetitions,1)
start_time <- Sys.time()
mult1=c(0.5)
n1=c(1000)
for (s in 1:1) {
  
  n=n1[s]
  for (l in 1:1) {
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
      
      
      ##creating the Delta_i functions
      h_t<-function(x_i_j,t){
        aux=1
        for(i in 1:(d))
        {
          aux=sqrt(2)*sin(t*(pi/2)*x_i_j[i])*aux
          
        }
        return(aux)
      }
      
      
      #defining the b_t_is
      b_i_ts=matrix(NaN,nrow=25,ncol=nm)
      #b_t_is=rnorm(btjs_val,0,1)
      for(i in 1:nm)
      {
        for (t in 1:25) {
          b_i_ts[t,i]=rnorm(1,0,20)
        }
      }
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
          y_i_js[aux-ms[i]+j]=f(X_s[aux-ms[i]+j,])+xi_i(X_s[aux-ms[i]+j,],i)+epsilon_s[aux-ms[i]+j]
        }
      }
      y_i_js
      Weights
      
      
      train_indices <- sample(1:length(y_i_js), round(length(y_i_js)*0.75), replace = FALSE)
      
      # subset the data into training and test sets
      X_s_training <- X_s[train_indices, ]
      y_i_js_training <- y_i_js[train_indices]
      Weights_training <- Weights[train_indices]
      X_s_test <- X_s[-train_indices, ]
      y_i_js_test <- y_i_js[-train_indices]
      Weights_test <- Weights[-train_indices]
      
      #######################################
      #######################################
      ####### Executing genlasso ############
      #######################################
      #######################################
      #######################################
      #######################################
      ###### Executing ######################
      ###CReating lambda regim for cv##########
      penalty=c(0.005,0.05,0.1,0.5,1,1.5,2)
      penalty_opt=0.5#cv_mars(X_s_training,y_i_js_training,folds_number,number_of_neighb,penalty)
      fit=mars(X_s_training, y_i_js_training,Weights=Weights_training, penalty=penalty[penalty_opt])
      yy=predict(fit,X_s_test,penalty=penalty[penalty_opt])
      result_rep_mars[rep,r]=mse_new(yy,X_s_test,Weights_test)
      
      cp=c(0.000001,0.00001,0.0001,0.001,0.01)
      cp_opt=0.001#cv_rpart(X_s_training,y_i_js_training,folds_number,number_of_neighb,cp)
      fit1= rpart(y_i_js_training~., data = data.frame(X_s_training, y_i_js_training),weights = Weights_training,
                  control =rpart.control( cp[cp_opt]))
      fit.pruned = prune(fit1, cp = cp[cp_opt])
      pred_y = predict(fit.pruned, data.frame(X_s_test))
      pred_y=as.matrix(pred_y)
      result_rep_cart[rep,r]=mse_new(pred_y,X_s_test,Weights_test)
      
      
      nodesize=c(1,3,5,8,10)
      nodesize_opt=1#cv_randomForest(X_s_training,y_i_js_training,folds_number,number_of_neighb,nodesize)
      fit2=randomForest(x=X_s_training,y=y_i_js_training,Weights = Weights_training,  ntree = 500,nodesize=nodesize[nodesize_opt])
      predy1=predict(fit2,newdata =X_s_test,type="response",ntree = 500,nodesize=nodesize[nodesize_opt])
      predy1=as.matrix(predy1)
      result_rep_RandFor[rep,r]=mse_new(predy1,X_s_test,Weights_test)
      print(rep)
    }}}
end_time <- Sys.time()
end_time - start_time
