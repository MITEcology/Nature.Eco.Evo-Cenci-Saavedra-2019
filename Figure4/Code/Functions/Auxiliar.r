### To make next step prediction
Testing <- function(C, C0, X){
  c0 = C0[nrow(C0), ]
  J = C[[length(C)]]
  return(c0 + J%*%X)
}
Add_to_TS <- function(TS, x){
  return(rbind(TS, x))
}
TakeLag <- function(X, species.to.lag, num.lag){
  tmp = matrix(0, nrow(X), num.lag)
  tmp[,1] = X[,species.to.lag]
  tmp[1, 1] = NA
  tmp[2:nrow(X), 1] = X[1:(nrow(X) - 1), species.to.lag]
  if(num.lag > 1){
    for(lag. in 2:num.lag){
      tmp[,lag.] = X[,species.to.lag]
      tmp[1, lag.] = NA
      tmp[2:nrow(X), lag.] = tmp[1:(nrow(tmp) - 1), lag.-1]
    }
  }
  tmp
}
make.lagged.ts <- function(X,sp.lag.selection ){
  ### X = time series
  ### sp.lag is a vector whose entry are the lags of each variable
  ### e.g., sp.lag = c(x.lag, y.lag, ..., u.lag)
  s = list()
  for(i in 1:length(sp.lag.selection)){
    Lag.sp = TakeLag(X, original.Embedding[i], sp.lag.selection[i])
    s[[i]] = cbind(X[,original.Embedding[i]], Lag.sp)
  }
  X = do.call(cbind,s)
  ### Remove the NA
  X = X[-c(1:max(sp.lag.selection)),]
  ### Save the position of the unlagged variables
  original.col = c()
  for(k in 1:length(sp.lag.selection)){
    if(k == 1){ original.col = c(original.col, 1)}else{
      num.lags = sum(unlist(lapply(1:(k-1), function(x,X) X[x], sp.lag.selection)))
      original.col = c(original.col, k + num.lags )
    }
  }
  return(list(time.series = X, original.variables = original.col))
}

take.coeff <- function(X, col.to.extract, original.emb){
  ### To use when prediction are made using lagged variables
  ### Take as input the sequence X of Jacobian along the attractor
  ### and the species to look at
  ### return a new sequence of Jacobian of the interaction among those species
  m = lapply(1:length(X$J), function(t, M, specie) M$J[[t]][specie,specie], 
                X, col.to.extract)
  for(i in 1:length(m)){
    colnames(m[[i]]) = rownames(m[[i]]) =original.emb
  }
  return(m)
}
naive.forecast <- function(last.point.training, test.set){
  #### Return the naive forecast, i.e., the test set is the last point of the training set
  naive.pred = matrix(0,nrow(test.set), ncol(test.set))
  for(j in 1:ncol(naive.pred)){
    naive.pred[,j] = rep(last.point.training[j], nrow(naive.pred))  
  }
  return(compute.rmse.test(naive.pred, test.set))
}
### Compute the rmse between two multivariate time series
compute.rmse.train <- function(X, Y){
  X = X[-1,]
  rmse = c()
  for(i in 1:ncol(X)){
    combine_xy = cbind(X[,i], Y[,i])
    rmse = c(rmse, sqrt(mean(unlist(lapply(1:nrow(combine_xy), 
                                           function(x, A) (A[x,1] - A[x,2])^2, combine_xy)))))
  }
  return(mean(rmse))
}
compute.rmse.test <- function(X, Y){
  rmse = c()
  for(i in 1:ncol(X)){
    combine_xy = cbind(X[,i], Y[,i])
    rmse = c(rmse, sqrt(mean(unlist(lapply(1:nrow(combine_xy), 
                                           function(x, A) (A[x,1] - A[x,2])^2, combine_xy)))))
  }
  return(mean(rmse))
}

ReadTimeSeries <- function(Nome){
  X = as.matrix(read.table(Nome))
  colnames(X) =  LETTERS[1:ncol(X)]
  return(X)
}
Standardizza <- function(X){
  ### This return y = (x-meanx)/stdx
  for(i in 1:ncol(X)){
    X[,i] = (X[,i]- mean(X[,i]))/sd(X[,i])
  }
  return(X)
}
Standardizza.test <- function(X, Y){
  ### X = test set
  ### Y = training set
  ### This return y = (x-meanY)/stdY
  for(i in 1:ncol(X)){
    X[,i] = (X[,i]- mean(Y[,i]))/sd(Y[,i])
  }
  return(X)
}
###########################
#### Here you compute the quality of the forecast as mean correlation coefficient
#### And we set to zero all those forecast that predict an extinction
MeanCorrelation <- function(TS, X){
  rho  = c()
  for(i in 1:ncol(X)){
   rho = c(rho, cor(TS[,i], X[,i]))
  }
  return(mean(rho))
}