Exponential.Kernel <- function(dst, theta){
  dbar <- mean(dst)
  krnl <- exp(-theta*dst/dbar)
  return(krnl)
}
Epanechnikov.Kernel <- function(dst, theta){
  tt <- dst/theta * ifelse(abs(dst/theta) > 1, 0, 1)
  krnl <- 0.75*(1 - tt^2)
  return(krnl)
}
TriCubic.Kernel <- function(dst,theta){
  tt <- dst/theta * ifelse(abs(dst/theta) > 1, 0, 1)
  krnl <- (1 - tt^3)^3
  return(krnl)
}
Matern.Kernel <- function(dst,theta){
	krnl = (1 + sqrt(3)*dst/theta)*exp(-(sqrt(3)*dst)/theta)
	return(krnl)
}
