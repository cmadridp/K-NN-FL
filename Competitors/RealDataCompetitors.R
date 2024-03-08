rm(list = ls())
#library(zoo)
#library(splines)
library(fda)
#library(foreach)
#library(doParallel)
library(ks)
#library(tidyverse)
#library(GGally)
#library(fields)
#library(ggplot2)
#library(dplyr)
#library(tidyr)
#library(faux)
library(mvtnorm)
#library(ggplot2)
#library(ggpubr)
#library(ggthemes)
#library(extrafont)
#library(datasets)
#library(ggpubr)
#library(gridExtra)  


library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
#library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(lubridate)
library(fields)

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
      mse_penalty[k,i]=mse_new(yy,y_i_js_training_test,Weights_training_test)
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
      mse_cp[k,i]=mse_new(pred_y,y_i_js_training_test,Weights_training_test)
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
      mse_nodesize[k,i]=mse_new(predy1,y_i_js_training_test,Weights_training_test)
    }
  }
  for(i in 1:(length(new_nodesize))){
    nodesize_min[i]=sum(mse_nodesize[,i])/folds_number
  }
  nodesize_opt=which(nodesize_min == min(nodesize_min))
  return(nodesize_opt)
}




mse_new <- function(y_predicted, y_i_js_test, Weights_tes) {
  diff <- y_predicted[,1] - y_i_js_test
  mse <- sum((diff^2) * Weights_tes)
  return(mse)
}




#The one I want
data=stack("sst.mon.mean.nc")
data1=nc_open("sst.mon.mean.nc")
lat=ncvar_get(data1,'lat')
lon=ncvar_get(data1,'lon')
time=ncvar_get(data1,'time')
time2=as.Date('1891-01-01')+time
sst=ncvar_get(data1,"sst")


m1=6 #June
i.year=1850
f.year=2019
t.years=f.year-i.year+1
iyear.pos=i.year-1850+1

m1.ipos=(iyear.pos-1)*12+m1 #m1 from i.year
m1.indices=sequence(t.years,from=m1.ipos,by=12)

m1.data=data[[m1.indices]]





#Region between europe and america in the pacific

image.plot(lon,rev(lat),sst[,length(lat):1,m1.indices[170]],xlab='Longitud',ylab='Latitude')
title("Sea Surface Temperature June 2019")
rect(xleft=300.5,xright = 339.5,ybottom=10.5,ytop=39.5, density=0, col = "Black",lwd=4) #Regions of interest corresponding to the Caribbean Sea 
rect(xleft=165.5,xright = 224.5,ybottom=35.5,ytop=54.5, density=0, col = "black",lwd=4) #Regions of interest corresponding to the Caribbean Sea 
rect(xleft=50.5,xright = 110.5,ybottom=-40.5,ytop=-20.5, density=0, col = "black",lwd=4) #Regions of interest corresponding to the Caribbean Sea 

#The region of interest is a 10x5 grid, i.e. n=50
lat_grid=seq(10.5,39.5,1)
lon_grid=seq(300.5,339.5,1)

total=length(lat_grid)*length(lon_grid)
matrix_lon_lat=matrix=matrix(NaN,total,2)
for(i in 0:(length(lat_grid)-1)){
  matrix_lon_lat[(i*length(lon_grid)+1):((i+1)*length(lon_grid)),]=lon_grid
  matrix_lon_lat[(i*length(lon_grid)+1):((i+1)*length(lon_grid)),2]=lat_grid[i+1]
}


aux_matrix_sst=matrix(NaN,total,170)
for(i in 1:170){
  for(j in 1:length(lat_grid)){
    for(k in 1:length(lon_grid)){
      if(is.na(sst[301:340,80:51,m1.indices[i]][k,j])){
        sst[301:340,80:51,m1.indices[i]][k,j]=0
      }
    }
  }
  
}


n=170
#the responses lies in dimension 3
m=total
#the observations lies in dimension 1
d=2
#we define the n_ts
nts<-matrix(NaN,nrow=n, ncol=1)
for(i in 1:n){
  nts[i]=m
}
nm_1=sum(nts)
nm_1


#creating X_s
X_s_b=matrix(,0,2)
for(aux in 1:170){
  X_s_b=rbind(X_s_b,matrix_lon_lat)
}
#creating y_t_is
y_i_js_b=matrix(NaN,nm_1,1)
for(i in 0:169){
  aux=as.vector(sst[301:340,80:51,m1.indices[i+1]])
  #print(aux)
  y_i_js_b[((i)*m+1):((i+1)*m)]=aux
}


#creating weights
Weights_b=matrix(NaN,nrow=nm_1,ncol=1)
aux=0
for(i in 1:n)
{
  aux=nts[i]+aux
  for(j in 1:(nts[i]))
  {
    Weights_b[aux-nts[i]+j]=1/(n*nts[i])
    
  }
}
Weights_b
write(y_i_js_b,file="testA.text")
write(X_s_b,file="AE_XA.text")
write(Weights_b,file="AE_WeightsA.text")



n <- 170
m <- total
num_splits <- 10
split_size <- m %/% num_splits

ms<-matrix(NaN,nrow=n, ncol=1)
for(i in 1:n){
  ms[i]=split_size
}
nm=sum(ms)
nm
result_rep_mars=matrix(NaN,num_splits,1)
result_rep_cart=matrix(NaN,num_splits,1)
result_rep_RandFor=matrix(NaN,num_splits,1)
for (s in 1:num_splits) {
  # get the current split indices
  start_index <- ((s-1)*split_size) + 1
  end_index <- s * split_size
  index_mask=rep(FALSE, nrow(X_s_b))
  # create index mask for the data within the current split
  for(j in 1:n){
    index_mask[(start_index+(j-1)*total):(end_index+(j-1)*total)] <- c(rep(TRUE, split_size))
  }
  # get the data for the current split
  X_s <- X_s_b[index_mask, ]
  y_i_js <- y_i_js_b[index_mask, ]
  Weights <- Weights_b[index_mask, ]
  
  library(caret)
  
  # set the seed for reproducibility
  #set.seed(123)
  
  # create the training and test indices
  train_indices <- createDataPartition(y_i_js, p = 0.75, list = FALSE)
  
  # subset the data into training and test sets
  X_s_training <- X_s[train_indices, ]
  y_i_js_training <- y_i_js[train_indices]
  Weights_training=Weights[train_indices]
  X_s_test <- X_s[-train_indices, ]
  y_i_js_test <- y_i_js[-train_indices]
  Weights_test=Weights[-train_indices]
  
  
  
  
  #######################################
  #######################################
  ####### Executing genlasso ############
  #######################################
  #######################################
  #######################################
  #######################################
  ###### Executing ######################
  ###CReating lambda regim for cv##########
  penalty=c(exp(-15),10,100,0.01,1)
  penalty_opt=1#cv_mars(X_s_training,y_i_js_training,folds_number,number_of_neighb,penalty)
  
  
  
  
  # fit the earth model using the training data
  
  
  
  
  # X_all=cbind(X_s_training,y_i_js_training)
  #df <- as.data.frame(X_all)
  #fit <- earth(y_i_js_training ~ V1+V2, data=df,penalty = penalty[penalty_opt])
  #yy <- predict(object=fit,newdata=X_s_test)
  
  
  fit=mars(X_s_training, y_i_js_training, penalty=penalty[penalty_opt])
  yy=predict(fit,X_s_test,penalty=penalty[penalty_opt])
  result_rep_mars[s]=mse_new(yy,y_i_js_test,Weights_test)
  
  cp=c(0.000001,0.00001,0.0001,0.001,0.01)
  cp_opt=1#cv_rpart(X_s_training,y_i_js_training,folds_number,number_of_neighb,cp)
  fit1= rpart(y_i_js_training~., data = data.frame(X_s_training, y_i_js_training),
              control =rpart.control( cp[cp_opt]))
  fit.pruned = prune(fit1, cp = cp[cp_opt])
  pred_y = predict(fit.pruned, data.frame(X_s_test))
  pred_y=as.matrix(pred_y)
  result_rep_cart[s]=mse_new(pred_y,y_i_js_test,Weights_test)
  
  
  nodesize=c(1,3,5,8,10)
  nodesize_opt=1#cv_randomForest(X_s_training,y_i_js_training,folds_number,number_of_neighb,nodesize)
  fit2=randomForest(x=X_s_training,y=y_i_js_training,  ntree = 500,nodesize=nodesize[nodesize_opt])
  predy1=predict(fit2,newdata =X_s_test,type="response",ntree = 500,nodesize=nodesize[nodesize_opt])
  predy1=as.matrix(predy1)
  result_rep_RandFor[s]=mse_new(predy1,y_i_js_test,Weights_test)
}
