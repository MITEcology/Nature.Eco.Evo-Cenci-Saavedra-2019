ELNET_fit_ <- function(time.series, targ_col, Embedding, theta, lambda,alp){
  Edim <- length(Embedding)
  coeff_names <- sapply(colnames(time.series),function(x) paste("d", targ_col, "d", x, sep = ""))
  block <- cbind(time.series[2:dim(time.series)[1],targ_col],time.series[1:(dim(time.series)[1]-1),])
  block <- as.data.frame(apply(block, 2, function(x) (x-mean(x))/sd(x)))
  
  lib <- 1:dim(block)[1]
  pred <- 1:dim(block)[1]
  
  coeff <- array(0,dim=c(length(pred),Edim + 1))
  colnames(coeff) <- c('c0', coeff_names)
  coeff <- as.data.frame(coeff)
  
  for (ipred in 1:length(pred)){
    libs = lib[-pred[ipred]]
    q <- matrix(as.numeric(block[pred[ipred],2:dim(block)[2]]),
                ncol=Edim, nrow=length(libs), byrow = T)
    distances <- sqrt(rowSums((block[libs,2:dim(block)[2]] - q)^2))
    ### Kernel
    Krnl = match.fun(Regression.Kernel)
    Ws = Krnl(distances, theta)
    ############ Fit function
    x = as.matrix(block[libs,2:dim(block)[2]])
    y = as.matrix(block[libs,1])
    x = x[seq_along(y), ]
    y = y[seq_along(y)]
    Ws = Ws[seq_along(y)]
    x = Ws * cbind(1, x)
    y = Ws * y
    fit <- enet(x, y, lambda = lambda, normalize = TRUE, intercept = FALSE)
    coeff[ipred,] <- predict(fit, s = alp, type="coefficients", mode="fraction")$coefficients 
  }
  return(coeff)
}

ELNET_fit <- cmpfun(ELNET_fit_)
