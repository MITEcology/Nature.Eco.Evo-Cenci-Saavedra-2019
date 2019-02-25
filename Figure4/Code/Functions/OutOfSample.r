out_of_sample_sequence<- function(cl, Jacobiano, th, lm, ts_training, num_points, ...){
  out_of_samp = c()
  coeff = list()
  ### Take the last point in the training set
  new_point = ts_training[nrow(ts_training), ]
  for(j in 1:num_points){
    ### Predict the first point in the training set and then allthe others
    new_point = Testing(Jacobiano$J, Jacobiano$c0, new_point)
    out_of_samp = rbind(out_of_samp, t(new_point))
    ts_training = Add_to_TS(ts_training, t(new_point))
    Jacobiano = update_Jacobian_TS(cl, ts_training, TargetList, Embedding, th, lm, RegressionType,alpha)
    coeff[[j]] = Jacobiano$J[[length(Jacobiano$J)]]
  }
  return(list(out_of_samp = out_of_samp, coeff = coeff))
}
####### The next functions serve to update the interaction without running cross validation
update <- function(cl, X, TargetList, Embedding, th, lm, RegressionType,alpha){
  J = c0 = list()
  n_ = 1
  FUN = match.fun(RegressionType)
  for(df in TargetList){
    ########## Now compute the optimum regularized coefficients
    J[[n_]]  = FUN(X, df, Embedding, th[n_], lm[n_],alpha)
    c0[[n_]] = J[[n_]]$c0
    J[[n_]] = J[[n_]][-1]
    n_ = n_ + 1
  }
  return(list(J = J, c0 = c0))
}
######################################################################
update_Jacobian_TS <- function(cl, X, TargetList, Embedding, th, lm, RegressionType,alpha){
  mine_output = update(cl, X, TargetList, Embedding, th, lm, RegressionType,alpha)
  mine_c0  = mine_output$c0
  mine_output = mine_output$J
  
  J = list()
  c0 = do.call(cbind, lapply(1:ncol(X), function(x, M) unlist(M[[x]]), mine_c0))
  colnames(c0) = sapply(TargetList,function(x) paste("c0_", x, sep = ""))
  for(k in 1:(nrow(X) - 1)){
    J[[k]] = do.call(rbind, lapply(1:ncol(X), function(x, M, i) unlist(M[[x]][i,]), mine_output, k))
    rownames(J[[k]]) = LETTERS[1:ncol(X)]
    colnames(J[[k]]) = LETTERS[1:ncol(X)]
    
  }
  return(list(J = J, c0 = c0))
}









