jacobiano_analitico <- function(nome, Inferred, num_species, origin, end){
  ### Inferred is the infer jacobian coefficients
  ### example jacobiano_analitico(BestCoefficients.elnet$J, 2, t.min, t.min+length.training-2)
  nomi = list()
  nomi$row = nomi$col = LETTERS[1:num_species]
  all.data = as.matrix(read.table(nome, header = F))
  ### Nota importante: e fondamentale che transponi la matrice per come li stampi da python
  J = lapply(origin:end, function(x, X) t(matrix(X[x,], num_species, 
                                                 num_species, dimnames = nomi)), all.data)
  result = Models.coefficient.correlation(Inferred, J, dfdx)
  tr.true = unlist(lapply(2:(length(J)), function(x, X) sum(diag(X[[x]])), J))
  tr.inferred = unlist(lapply(1:length(Inferred), function(x, X) sum(diag(X[[x]])), Inferred))
  tr.inferred = tr.inferred[!is.na(tr.inferred)]
  trace.corr = cor(scale(tr.true), scale(tr.inferred))
  return(list(correlation.matrix = result$mat.cor, rmse.matrix = result$mat.rmse, 
              analytical.jacobian = J, trace.correlation = trace.corr,
              true.trace = tr.true, inferred.trace = tr.inferred))
}

Models.coefficient.correlation <- function(X, Y, dfdx){
  ### input X,Y = Best series of jacobian coefficients
  ###       dfdx = a dataframe with all the combination of names of species (Jacobian entry)
  ### output: A matrix of which each entry tells me the correlation coefficient between 
  ###         matrix X and matrix Y
  CorMat = rmse.mat = c()
  for(k in 1:nrow(dfdx)){
    X_ij = unlist(lapply(1:(length(X)), function(i, M) M[[i]][dfdx[k,1],dfdx[k,2]], X))
    Y_ij = unlist(lapply(1:(length(Y)-1), function(i, M) M[[i]][dfdx[k,1],dfdx[k,2]], Y))
    CorMat = c(CorMat, cor(X_ij, Y_ij))
    rmse.mat = c(rmse.mat, sqrt(mean((X_ij - Y_ij)^2)))
  }
  return(list( mat.cor = matrix(CorMat, nrow(X[[1]]), nrow(X[[1]])),
               mat.rmse = matrix(rmse.mat,nrow(X[[1]]), nrow(X[[1]]))))
}

