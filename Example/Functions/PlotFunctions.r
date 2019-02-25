############ For the plotting of the time series
Put_train_test_together <- function(d.tr, d.ts, pr.tr, pr.ts, species){
  title = paste('Species', species)
  
  plot(c(d.tr[,species],d.ts[,species]), pch = 20, ylim = range(c(d.tr[,species]-1, pr.tr[,species]+1)), 
       main = title, ylab = expression('x'[species]))
  lines(c(pr.tr[,species], pr.ts[,species]), col = 'red', lwd = 4, lty = 1)
  abline(v=(nrow(d.tr)+1), col="blue", lty = 2, lwd = 1.5)
  text(nrow(d.tr)/2, 2.5,'Training set')
  text(nrow(d.tr) + 8, 2.5,'Test set')
}