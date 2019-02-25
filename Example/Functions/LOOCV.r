########################### Cross Validation ##################################
LOO.CV <- function(cl, X, TargetList, Embedding, grid, RegressionType,alpha, Regression.Kernel){
  mine_output = Jacobian_(cl, X, TargetList, Embedding, grid, RegressionType,alpha, Regression.Kernel)
  validation.error = mine_output$min_validation_error
  theta_opt = mine_output$th
  lambda_opt = mine_output$lm
  mine_c0  = mine_output$c0
  mine_output = mine_output$J
  
  J_ = list()
  C0_ = do.call(cbind, lapply(1:ncol(X), function(x, M) unlist(M[[x]]), mine_c0))
  colnames(C0_) = sapply(TargetList,function(x) paste("c0_", x, sep = ""))
  for(k in 1:(nrow(X) - 1)){
    J_[[k]] = do.call(rbind, lapply(1:ncol(X), function(x, M, i) unlist(M[[x]][i,]), mine_output, k))
    rownames(J_[[k]]) = LETTERS[1:ncol(X)]
    colnames(J_[[k]]) = LETTERS[1:ncol(X)]
    
  }
  BestCoefficients = list()
  BestCoefficients$J = J_
  BestCoefficients$c0 = C0_
  BestParameters = list()
  BestParameters$BestTH = theta_opt
  BestParameters$BestLM = lambda_opt
  return(list(BestCoefficients = BestCoefficients, BestParameters = BestParameters, validation.error = validation.error))
}

##########
LOOCrossValidation <- function(cl, data.df, targ.sp, Embedding, grid, RegressionType,alpha, Regression.Kernel){
  S_target <- parLapply(cl, 1:nrow(grid), CV, data.df, targ.sp, Embedding, grid, RegressionType,alpha, Regression.Kernel)
  error.mat = cbind(unlist(S_target)[attr(unlist(S_target),"names") == "bndwth"],
                    unlist(S_target)[attr(unlist(S_target),"names") == "lmb"],
                    unlist(S_target)[attr(unlist(S_target),"names") == "mse"])
  rownames(error.mat) = c()
  error.mat[,3] = round(error.mat[,3],4)
  error.mat = error.mat[order(error.mat[,3]), ]
  idx = which(error.mat[,3] == min(error.mat[,3]))
  idx.th = which(grid[idx,1] == min(grid[idx,1]))
  idx = idx[min(idx.th)]
  return(list(BestTH = error.mat[idx,1],
              BestLM = error.mat[idx,2],
              min.val.err = error.mat[idx,3],
              val.err = error.mat))
}

#################################################
CV <- function(i, time_series_training, target_species, Embedding, TL, RegressionType,alpha, Regression.Kernel){
  FUN = match.fun(RegressionType)
  #### Fit the model
  coefficients = FUN(time_series_training, target_species, Embedding, TL[i,1], TL[i,2],alpha, Regression.Kernel)
  #### Take the variables in the time series that belong to the embedding
  time_series_training = time_series_training[,Embedding]
  #### Standardize the time series
  time_series_training = Standardizza(time_series_training)
  #### Take the forecast
  Data = Forecast_SMap(coefficients, time_series_training)
  time_series_training = time_series_training[-1,]
  Data = Data[-length(Data)]
  D = cbind(time_series_training[,target_species], Data)
  #### Compute the mean square error
  MSE = mean(unlist(lapply(1:nrow(D), function(x, A) (A[x,1] - A[x,2])^2, D)))
  #return(MSE)
  return(list(mse = MSE, bndwth = TL[i,1], lmb = TL[i,2]))
}

######################################################
Forecast_SMap <- function(C, X){
  n_coeff = ncol(X)
  predizione = c()
  for(k in 1:nrow(X)){
    s = unlist(lapply(1:n_coeff, function(x, J, i) return(J[[x + 1]][i]), C, k))
    predizione = c(predizione, C[[1]][k] + s%*%X[k,])
  }
  return(predizione)
}

Jacobian_ <- function(cl, X, TargetList, Embedding, grid, RegressionType,alpha, Regression.Kernel){
  J = c0 = list()
  th = lm = val.error = c()
  n_ = 1
  FUN = match.fun(RegressionType)
  for(df in TargetList){
    RegularizedParameters <- LOOCrossValidation(cl, X, df, Embedding, grid, RegressionType,alpha, Regression.Kernel)
    ########## Now compute the optimum regularized coefficients
    J[[n_]]  = FUN(X, df, Embedding, RegularizedParameters$BestTH, RegularizedParameters$BestLM,alpha, Regression.Kernel)
    th = c(th, RegularizedParameters$BestTH)
    lm = c(lm, RegularizedParameters$BestLM)
    c0[[n_]] = J[[n_]]$c0
    J[[n_]] = J[[n_]][-1]
    n_ = n_ + 1
    val.error = c(val.error, RegularizedParameters$min.val.err)
  }
  min_validation_error = mean(val.error)
  return(list(J = J, c0 = c0, th = th, lm = lm, min_validation_error = min_validation_error))
}
