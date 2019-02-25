suppressMessages(library(Hmisc))
options(warn=-1)
#######################################################################################
######### Compute the structural identifiability of ecological interactions  ##########
#######################################################################################
wd.name = getwd()
ModelName = 'NewZeland'
file.name = paste(wd.name, '/Results/result_from_', ModelName, '.RData', sep =  '')
data.name = paste(wd.name, '/Data/data_from_', ModelName, '.RData', sep = '')
load(data.name)
rm(data.name,wd.name)
##### round the error to two digits
epsilon.round.rmse = 2
### subset those models with error at maximum 10% bigger than the minimum
idx.rmse.test = which(round(test.error.rmse,epsilon.round.rmse) < 
                        1.05*min(round(test.error.rmse,epsilon.round.rmse))) 
### subset those models with explained variance at minimum 90% of the maximum explained variance
idx.r2.train = which(round(Explained.Variance.train,epsilon.round.rmse) > 
                        0.9*max(round(Explained.Variance.train[idx.rmse.test],epsilon.round.rmse)))
### Find the intersection between the best models in the training and test set
idx.subset = intersect(idx.r2.train, idx.rmse.test)
####################################################################################################
species.names = colnames(jac.coefficients[[1]][[1]])
idx.to.use.for.subset = idx.subset
#############################################
##### Coefficients of the ensemble method with their standard error
##### ensemble.coefficients will be the time series of Jacobian from the ensemble method
##### ensemble.coefficients.se is the time series of their standard errors
##### ensemble.coefficients.cv is the the series of coefficient of variations
ensemble.coefficients = ensemble.coefficients.se =  ensemble.coefficients.cv = list()
##### The weights are like probability of a model and must sum up to one
metodi = c('R2', 'uniform', 'all.models')
methods.pesi = metodi[1]
if(methods.pesi == 'R2'){
    pesi = Explained.Variance.train[idx.to.use.for.subset]/sum(Explained.Variance.train[idx.to.use.for.subset])
    cat('weights chosen from', methods.pesi, '\n')
} else if(methods.pesi == 'uniform'){
    pesi = rep(1./length(idx.to.use.for.subset), length(idx.to.use.for.subset))
    cat('weights chosen from', methods.pesi, '\n')
} else if(methods.pesi == 'all.models'){
    idx.to.use.for.subset = 1:length(test.error.rmse)
    pesi = Explained.Variance.train[idx.to.use.for.subset]/sum(Explained.Variance.train[idx.to.use.for.subset])
    cat('weights chosen from', methods.pesi, '\n')
}
### Compute ensemble coefficients, the standard errors and the coefficient of variations
for(k in 1:length(jac.coefficients[[1]])){
  ensemble.coefficients[[k]] = matrix(NA,num.species,num.species)
  ensemble.coefficients.se[[k]] = matrix(NA,num.species,num.species)
  ensemble.coefficients.cv[[k]] = matrix(NA,num.species,num.species)    
  for(i in 1:num.species){
    for(j in 1:num.species){
      boosted.coeff = unlist(lapply(idx.to.use.for.subset, function(x,X,l,m,n) X[[x]][[n]][l,m], 
                                    jac.coefficients,i,j,k))
  #### Weigthed mean over models
      ensemble.coefficients[[k]][i,j] = wtd.mean(boosted.coeff, pesi)
	#### Weighted standard error
      ensemble.coefficients.se[[k]][i,j] = 
      1.96*sqrt(wtd.var(boosted.coeff, weights = pesi, na.rm = T,normwt=T))/sqrt(length(idx.to.use.for.subset))
	#### Weighted coefficient of variation
	    ensemble.coefficients.cv[[k]][i,j] = 
        sqrt(wtd.var(boosted.coeff, weights = pesi, na.rm = T,normwt=T))/abs(wtd.mean(boosted.coeff, weights = pesi))
    }
  }
}
subsetted.divergence = list()
count = 1
for(dvg in idx.to.use.for.subset){
  subsetted.divergence[[count]] = divergence[[dvg]]
  count = count + 1
}
traccia = traccia.error = rep(0,length(subsetted.divergence[[1]]))
subsetted.divergence = sapply(subsetted.divergence, unlist)
for(s in 1:nrow(subsetted.divergence)){
      traccia[s] = wtd.mean(subsetted.divergence[s,], pesi)
      diagonal.elements = diag(ensemble.coefficients.se[[s]])
      ### Compute the error by propagating the error on the element of the diagonal
      traccia.error[s] = sqrt(sum(unlist(lapply(1:length(diagonal.elements), function(x,X) X[x]^2, diagonal.elements))))
}
###########################
save(species.names, 
     ensemble.coefficients, ensemble.coefficients.se, ensemble.coefficients.cv,
     traccia, traccia.error,
     file = file.name)

