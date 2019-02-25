Reconstructing <- function(J, C0, X){
  return(C0 + matrix(unlist(J),length(X), length(X)) %*% X)
}
###
RootMeanSquareError <- function(X, Y){
rmse <- unlist(lapply(1:length(X), function(x, X, Y) (X[x] - Y[x])^2, X, Y))
return(sqrt(mean(rmse)))
}
###
ReconstructionOfTrainingSet <- function(TimeSeries, Jac){
  prd_ = c()
  TimeSeries = Standardizza(TimeSeries)
  for(i in 1:(nrow(TimeSeries) - 1)){
    prd_ <- rbind(prd_, t(Reconstructing(Jac$J[i], Jac$c0[i,], TimeSeries[i,])))
  }
  colnames(prd_) = colnames(TimeSeries)
  return(prd_)
}
#####################################################################
#####################################################################
ComputeTrainingError <- function(TimeSeries, Jac, real.var){
  ############ First I check how the reconstructed time series looks like
  rmse = rho = c()
  TimeSeries = Standardizza(TimeSeries)
  prd_ = ReconstructionOfTrainingSet(TimeSeries, Jac)
  colnames(prd_) = colnames(TimeSeries)
  ########################################################################
  for(i in real.var){
    #### Important: you need to remove the first point because you 
    #### are making in sample predictions
    observed = TimeSeries[,i]
    observed = observed[-1]
    rmse = c(rmse, RootMeanSquareError(prd[,i], observed))
    rho = c(rho, cor(prd_[,i], observed))
  }
  return(list(rmse = median(rmse), rho = median(rho)))
}



