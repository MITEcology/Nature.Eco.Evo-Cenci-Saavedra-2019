############## Some plot function
###### To plot the inference of the coefficients
PlotCorrelationBetweenAlgorithms <- function(x,xlabel,ylabel){
  
  column_name <- sapply(Embedding,function(x) paste("d", x, sep = ""))
  xval <- formatC(x, format="f", digits= 3)
  pal <- colorRampPalette(c(rgb(0.96,0.96,1), rgb(0.05,0.1,0.9)), space = "rgb")
  hm <- heatmap.2(x, Rowv=FALSE, Colv=FALSE, dendrogram="none", 
                  col=pal, tracecol="#303030", trace="none", 
                  xlab = xlabel, ylab = ylabel, labRow = column_name, labCol = column_name,
                  cellnote=xval, notecol="black", notecex=0.7,  keysize = 0.2, key = FALSE)
  
}

PlotInferenceMatrix <- function(x, MainTitle){
  
  column_name <- sapply(Embedding,function(x) paste("d", x, sep = ""))
  xval <- formatC(x, format="f", digits= 3)
  pal <- colorRampPalette(c(rgb(0.96,0.96,1), rgb(0.05,0.1,0.9)), space = "rgb")
  
  hm <- heatmap.2(x, Rowv=FALSE, Colv=FALSE, dendrogram="none", # main = MainTitle,
                  col=pal, tracecol="#303030", trace="none", 
                  labRow = column_name, labCol = column_name,
                  cellnote=xval, notecol="black", notecex=0.7,  keysize = 0.2, key = FALSE)
  
}

############ For the plotting of the time series
Plot_out <-function(Xtest, Y){
  par(mfrow=c(1,ncol(Xtest)))
  for(k in 1:ncol(Xtest)){
    title = paste('Species', k)
    plot(Xtest[,k], pch = 20, ylim = range(c(Xtest[,k], Y[,k])), 
         main = title, ylab = expression('x'[k]))
    lines(Y[,k], col = 'red', lwd = 4)
    lines(Xtest[,k], lwd = 4)
  }
  par(mfrow=c(1,1))
}

Plot_AllSpeciesTraining <-function(Xtr, Yrec){
  par(mfrow=c(1,ncol(Xtr)))
  for(k in 1:ncol(Xtr)){
    title = paste('Species', k)
    plot(Xtr[,k], pch = 20, ylim = range(c(Xtr[,k], Yrec[,k])), 
         main = title, ylab = 'population size', xlab = 'Time')
    lines(Yrec[,k], col = 'red', lwd = 4, lty = 2)
    lines(Xtr[,k], lwd = 4)
  }
  
  par(mfrow=c(1,1))
}
Plot_SingleSpeciesTraining <- function(Xtr, Yrec, species){
  title = paste('Species', species)
  plot(Xtr[,species], pch = 20, ylim = range(c(Xtr[,species], Yrec[,species])), 
       main = title, ylab = 'population size', xlab = 'Time')
  lines(Yrec[,species], col = 'red', lwd = 4, lty = 2)
  lines(Xtr[,species], lwd = 4)
}