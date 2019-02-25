#############
structural.identifiability <- function(file.name, plt = TRUE){
  load(file.name)
  if(file.exists(file.name)){
    load(file.name)
    coeff.var.dist = list()
    tmp = 1
    num.species = ncol(ensemble.coefficients[[1]])
    if(plt == TRUE){par(mfrow=c(num.species,num.species))}
    for(l in 1:num.species){
      for(m in 1:num.species){
        x = l
        y = m
        coeff = confidence = coef.var = rep(0,length(ensemble.coefficients))
        for(i in 1:length(ensemble.coefficients)){
          coeff[i] = ensemble.coefficients[[i]][x,y]
          confidence[i] = ensemble.coefficients.se[[i]][x,y]
          coef.var[i] = ensemble.coefficients.cv[[i]][x,y]
        }
        element.name = paste('d', species.names[x], 'd', 
              species.names[y], sep = '')
        coeff.var.dist[[tmp]] = coef.var
        names(coeff.var.dist[[tmp]]) = element.name
        tmp = tmp + 1
        if(plt == TRUE){
          plot(coeff, type = 'l', lwd = 1, col = 'red', xlab = 'Time', 
               main = element.name)
          lines(coeff+(confidence), lwd  =1)
          lines(coeff-(confidence), lwd  =1)
        }
      }
    }
    if(plt == TRUE){par(mfrow=c(1,1))}
  }
  return(coeff.var.dist)
}
