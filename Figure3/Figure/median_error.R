rm(list = ls())
mdn = std = c()
conf95 <- function(X){
  wx = wilcox.test(X, alternative = "two.sided", conf.int = T, conf.level = 0.95)
  return((wx$conf.int[2] - wx$conf.int[1])/2)
}

a = as.matrix(read.table('lv.txt'))[,4]
mdn = c(mdn, median(na.omit(a)))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(na.omit(a), 1000, theta))
std = c(std, conf95(a))
################################################
a = as.matrix(read.table('fc.txt'))[,4]
mdn = c(mdn, median(a))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(a, 1000, theta))
std = c(std, conf95(a))
################################################
a = as.matrix(read.table('cr.txt'))[,4]
mdn = c(mdn, median(a))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(a, 1000, theta))
std = c(std, conf95(a))
################################################
a = as.matrix(read.table('ml.txt'))[,4]
mdn = c(mdn, median(a))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(a, 1000, theta))
std = c(std, conf95(a))
###############################################
a = as.matrix(read.table('hs.txt'))[,4]
mdn = c(mdn, median(a))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(a, 1000, theta))
std = c(std, conf95(a))
################################################
a = as.matrix(read.table('k.txt'))[,4]
mdn = c(mdn, median(a))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(a, 1000, theta))
std = c(std, conf95(a))
################################################
a = as.matrix(read.table('a.txt'))[,4]
mdn = c(mdn, median(a))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(a, 1000, theta))
std = c(std, conf95(a))
################################################
a = as.matrix(read.table('d.txt'))[,4]
mdn = c(mdn, median(a))
theta <- function(x){median(x)} 
#std = c(std, BootstrapError(a, 1000, theta))
std = c(std, conf95(a))
################################################
to.print = cbind(1:8, mdn, std)
write.table(to.print, file = 'median_plus_error.txt', row.names = F, col.names = F)
