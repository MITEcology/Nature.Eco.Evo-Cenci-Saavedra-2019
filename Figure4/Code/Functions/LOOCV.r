########################### Cross Validation ##################################
BestModelLOOCV <- function(cl, X, TargetList, Embedding, grid, RegressionType,alpha){
  mine_output = Jacobian_(cl, X, TargetList, Embedding, grid, RegressionType,alpha)
  theta_opt = mine_output$th
  lambda_opt = mine_output$lm
  mine_c0  = mine_output$c0
  mine_output = mine_output$J
  
  J_ = list()
  C0_ = do.call(cbind, lapply(1:ncol(X), function(x, M) unlist(M[[x]]), mine_c0))
  colnames(C0_) = sapply(TargetList,function(x) paste("c0_", x, sep = ""))
  for(k in 1:(nrow(X) - 1)){
    J_[[k]] = do.call(rbind, lapply(1:ncol(X), function(x, M, i) unlist(M[[x]][i,]), mine_output, k))
    rownames(J_[[k]]) = Embedding
    colnames(J_[[k]]) = Embedding
    
  }
  BestCoefficients = list()
  BestCoefficients$J = J_
  BestCoefficients$c0 = C0_
  BestParameters = list()
  BestParameters$BestTH = theta_opt
  BestParameters$BestLM = lambda_opt
  return(list(BestCoefficients = BestCoefficients, BestParameters = BestParameters))
}

##########
LOOCrossValidation <- function(cl, data.df, targ.sp, Embedding, grid, RegressionType,alpha){
  #### It is a bit long but at least you are sure that paralelizzation does not mess things up
  S_target <- parLapply(cl, 1:nrow(grid), CV, data.df, targ.sp, Embedding, grid, RegressionType,alpha)
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
              min.var.err = error.mat[idx,3],
              val.err = error.mat))
}
LOOCV.not.parallel <- function(cl, data.df, targ.sp, Embedding, grid, RegressionType){
  S_target <- unlist(lapply(1:nrow(grid), CV, data.df, targ.sp, Embedding, grid, RegressionType,alpha))
  S_target = round(S_target,4)
  idx = which(S_target == min(S_target))
  idx.th = which(grid[idx,1] == min(grid[idx,1]))
  idx = idx[min(idx.th)]
  return(list(BestTH = grid[idx,1], 
              BestLM = grid[idx,2], 
              min.val.error = min(S_target),
              ValidationError = matrix(S_target, attributes(grid)$out.attrs$dim[1], 
                           attributes(grid)$out.attrs$dim[2])))
}

#################################################
CV <- function(i, time_series_training, target_species, Embedding, TL, RegressionType,alpha){
  FUN = match.fun(RegressionType)
  #### Fit the model
  coefficients = FUN(time_series_training, target_species, Embedding, TL[i,1], TL[i,2],alpha)
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

Jacobian_ <- function(cl, X, TargetList, Embedding, grid, RegressionType,alpha){
  J = c0 = list()
  th = lm = c()
  n_ = 1
  FUN = match.fun(RegressionType)
  for(trg in TargetList){
    RegularizedParameters <- LOOCrossValidation(cl, X, trg, Embedding, grid, RegressionType,alpha)
    ########## Now compute the optimum regularized coefficients
    J[[n_]]  = FUN(X, trg, Embedding, RegularizedParameters$BestTH, RegularizedParameters$BestLM,alpha)
    th = c(th, RegularizedParameters$BestTH)
    lm = c(lm, RegularizedParameters$BestLM)
    c0[[n_]] = J[[n_]]$c0
    J[[n_]] = J[[n_]][-1]
    n_ = n_ + 1
  }
  return(list(J = J, c0 = c0, th = th, lm = lm))
}
