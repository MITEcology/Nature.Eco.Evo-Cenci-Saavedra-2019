### To make next step prediction
Testing <- function(C, C0, X){
  c0 = C0[nrow(C0), ]
  J = C[[length(C)]]
  return(c0 + J%*%X)
}
Add_to_TS <- function(TS, x){
  return(rbind(TS, x))
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
  X = as.matrix(read.table(Nome, header = F))
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
MeanCorrelation <- function(TS, X){
  rho  = rmse = c()
  for(i in 1:ncol(X)){
    rho = c(rho, cor(TS[,i], X[,i]))
    rmse = c(rmse, sqrt(mean((TS[,i] - X[,i])^2)))
  }
  return(list( correlation = mean(rho), rmse = mean(rmse)))
}
