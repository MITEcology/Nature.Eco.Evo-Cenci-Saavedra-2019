rm(list=ls())
source('functions.r')
stampa = FALSE
############################################# Data explanation
#################### path to data and results
wd.name = getwd()
ModelName = 'NewZeland'
data.data = paste('Data/data_from_', ModelName, '.RData', sep = '')
results.data = paste('Results/result_from_', ModelName, '.RData', sep = '')
fil.name = paste('structural_identifiability_', ModelName, '.txt', sep = '')
################################# Now look at the confidence intervals
dist.of.coef.var =structural.identifiability(results.data, plt = F)
nomi.interactions = unlist(lapply(1:length(dist.of.coef.var), function(x,X)
  names(X[[x]][1]), dist.of.coef.var))
################################# Now look at the coefficients of variation
##### Structural identifiability of coefficients
dist.of.coef.var = do.call(cbind, dist.of.coef.var)
boxplot(dist.of.coef.var, ylim = c(0,6.))
colnames(dist.of.coef.var) = nomi.interactions
rownames(dist.of.coef.var) = NULL
if(isTRUE(stampa)){
  for(i in 1:ncol(dist.of.coef.var)){
    to.print = na.omit(dist.of.coef.var[,i])
    fil.name = paste('boxplot_', i, '.txt', sep = '')
    write.table(to.print, file =fil.name , row.names = F, col.names = F, na = ' ')
  }
}
##### Now plot the trace with its error
load(results.data)
plot(traccia, type = 'l', lwd = 3, ylim = range(c(traccia + traccia.error, traccia - traccia.error)),
     xlab = 'Time', ylab = 'Divergence of vector field')
lines(traccia + traccia.error, lwd = .5, col = 'red')
lines(traccia - traccia.error, lwd = .5, col = 'red')
#####
library(lubridate)
empirical.data = as.matrix(read.table('NewZeland.txt', header = T))
time.interval = empirical.data[1:248,1]
mesi = month(as.POSIXlt(time.interval, format = "%d/%m/%Y"))
mt= cbind(mesi, traccia, traccia.error)
tmp = tmp.err = c()
for(i in 1:12){
  ### Subset months
  s = mt[mt[,1] == i, ]
  ### take the mean of the trace
  tmp = c(tmp, mean(s[,2]))
  ### prepagate the error
  tmp.err = c(tmp.err ,sqrt(sum(s[,3]^2))/length(s))
}
plot(tmp, type = 'h', xlab = 'Month', ylab = 'divergence of vector field', xaxt = "n")
axis(1, at=1:12, labels=c('Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov', 'Dec'))

if(isTRUE(stampa)){
  traccia.con.errore = cbind(1:length(traccia), traccia, traccia.error)
  colnames(traccia.con.errore) = c('Time', 'Trace', 'error')
  write.table(traccia.con.errore, file = '../../Figure/time_series_of_divergence.txt', row.names = F, na = ' ')
  ####
  data.to.print = cbind(1:12, tmp, tmp.err)
  colnames(data.to.print) = c('x', 'y', 'error')
  write.table(data.to.print, file = '../../Figure/Monthly_average_divergence.txt', row.names = F, na = ' ')
}