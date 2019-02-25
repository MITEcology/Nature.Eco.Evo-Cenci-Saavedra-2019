Reconstructing <- function(J, C0, X){
  return(C0 + matrix(unlist(J),length(X), length(X)) %*% X)
}
###
MeanSquareError <- function(X, Y){
mse <- unlist(lapply(1:length(X), function(x, X, Y) (X[x] - Y[x])^2, X, Y))
return(mean(mse))
}
###
ReconstructionOfTrainingSet <- function(TimeSeries, Jac){
  prd = c()
  TimeSeries = Standardizza(TimeSeries)
  for(i in 1:(nrow(TimeSeries) - 1)){
    prd <- rbind(prd, t(Reconstructing(Jac$J[i], Jac$c0[i,], TimeSeries[i,])))
  }
  return(prd)
}
#####################################################################
#####################################################################
ComputeTrainingError <- function(TimeSeries, Jac){
  ############ First I check how the reconstructed time series looks like
  mse = rho = c()
  TimeSeries = Standardizza(TimeSeries)
  prd = ReconstructionOfTrainingSet(TimeSeries, Jac)
  ########################################################################
  for(i in 1:ncol(TimeSeries)){
    observed = TimeSeries[,i]
    #### Important: you need to remove the first point because you 
    #### are making in sample predictions
    observed = observed[-1]
    mse = c(mse, MeanSquareError(prd[,i], observed))
    rho = c(rho, cor(prd[,i], observed))
  }
  return(list(MSE = median(mse), rho = median(rho)))
}



