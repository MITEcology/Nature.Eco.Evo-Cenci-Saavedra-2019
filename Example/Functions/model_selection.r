rm(list = ls())
suppressMessages(library(Matrix))
suppressMessages(library(quantreg))
suppressMessages(library(parallel))
suppressMessages(library(compiler))
suppressMessages(library(lars))
suppressMessages(library(elasticnet))
suppressMessages(library(caret))
options(warn=-1)
#####################################################################################
source('Auxiliar.r')
source('elastic_net_fit.r')
source('LOOCV.r')
source('KernelFunctions.r')
source('OutOfSample.r')
source('Interactions.R')
source('TrainingError.r')
################# Regularized S-map for inference and forecasting ###################
##### Few notes:
##### 1) Here we use elastic net regularization and three different kernel functions
##### 2) This is the most complete version of the code as it cross-validate 
#####    every possible variable:
#####     a) kernel function
#####     b) ratio of L1 to L2 norm
#####     c) kernel bandwidth
#####     d) lambda
##### For the analysis of empirical data with multivariate time series one should also
##### cross validate features and time lags
##### Generally, the optimum kernel function and ratio of L1 to L2 norm (as well as features 
##### and time lags in empirical data) are consistent across the time series.
##### So, if one were to repeat the inference and forecasting across multiple 
##### chunks of the time series. then a) and b) can be selected on a small training 
##### set and then kept fixed for the following analysis
##### On the other hand it is important to always cross validate c) and d)
##### -----------------------
##### Generally we suggest to run 'model_selection.r' (<--- very lengthy and slow)
##### Then get the optimum a) and b) and use them to run 'main.r' which is much faster
#####################################################################################
ShowPlot = F
options.models = c('inputFiles/PredatorPrey.txt','inputFiles/RPS.txt','inputFiles/Deterministic_TimeSeries.txt')
options.jacobian = c('inputFiles/jacobian_predator.txt','inputFiles/jacobian_rps.txt','inputFiles/jacobian_chaos.txt')
Kernel.Options = c('Exponential.Kernel', 'Epanechnikov.Kernel', 'TriCubic.Kernel')
choice = 1
cat('Running model selection for:', options.models[choice])
FileName = options.models[choice]
JacobianName = options.jacobian[choice]
###################################
logspace <- function(d1, d2, n) exp(log(10)*seq(d1, d2, length.out=n)) 
std_err <- function(x) sd(x)/sqrt(length(x))
####################################
   
##### length of training and test set
length.training = 100
length.testing = 10
### Read Time series
d = ReadTimeSeries(FileName)
### Use all the species
Embedding = LETTERS[1:ncol(d)]
TargetList = Embedding
######################
dfdx = expand.grid(TargetList, TargetList)
######################
d = d[, Embedding]
t.min = 1
####
ts.train.preserved = PreserveTimeSeriesForJacobian(d, length.training,
                                                   length.testing, t.min)$Xtr
#### Random Chunk of length:
length_of_interval = length.training + length.testing
t.max = t.min + length_of_interval - 1
interval = t.min:t.max
interval_training = 1:length.training
interval_testing = (length.training + 1):length_of_interval
#### Subset the chunk
d_intact = d[interval,]
d.training = d_intact[interval_training, ]
d.training = Standardizza(d.training)
#### Make noise if you want
ObservationalNoise = 0
if(ObservationalNoise == 1){
  d.training = d.training + d.training*matrix(rnorm(length(d.training), 0,0.1), 
                                              nrow(d.training), ncol(d.training))   
  cat('Now the training set is noisy\n')
}
#### You now need to standardize the test set using mean and sd of the training set
#### Notice that the test set is declared here but not used in the fit
d.testing = d_intact[interval_testing, ]
d.testing = Standardizza.test(d.testing,ts.train.preserved)
####
RegressionType = 'ELNET_fit'
ratio = c(0.85, 0.9, 0.95, 1.)
choice = c(1,2,3)
other.parameters = expand.grid(ratio,choice)
val_err = c()
#### Now run the fit in parallel
Lavoratori = detectCores() - 1
cl <- makeCluster(Lavoratori, type = "FORK")
for(k in 1:nrow(other.parameters)){
  Regression.Kernel = Kernel.Options[other.parameters[k,2]]
  lambda = logspace(-3,0,15)                       
  if(Regression.Kernel == 'Exponential.Kernel'){
    tht = logspace(-2,1.2,15)         
  }else{
    tht = logspace(-2,1.2,15)         
  }
  parameters_on_grid = expand.grid(tht, lambda)  
  BestModelElnet = BestModelLOOCV(cl, d.training, TargetList, Embedding, parameters_on_grid, 
                                  RegressionType,other.parameters[k,1],Regression.Kernel)
  val_err = c(val_err, BestModelElnet$validation.error)
}
idx = which(val_err == min(val_err))
Regression.Kernel = Kernel.Options[other.parameters[idx,1]]
lambda = logspace(-3,0,15)                       
if(Regression.Kernel == 'Exponential.Kernel'){
  tht = logspace(-2,1.2,15)         
}else{
  tht = logspace(-2,1.2,15)         
}
BestModelElnet = BestModelLOOCV(cl, d.training, TargetList, Embedding, parameters_on_grid, 
                                RegressionType,other.parameters[idx,1],Regression.Kernel)
stopCluster(cl)
#### Coefficient of the Jacobian and intercept
BestCoefficients.elnet = BestModelElnet$BestCoefficients
#### Bandwith of the kernel and regularization parameter
BestParameters.elnet = BestModelElnet$BestParameters
### Out-of-sample forecast
out.of.samp = out_of_sample_sequence(BestCoefficients.elnet, BestParameters.elnet$BestTH,
       	                                   BestParameters.elnet$BestLM, other.parameters[idx,1],
               	                           d.training, length.testing, Regression.Kernel)
prd = out.of.samp$out_of_samp
###############################
###### Now check the in-sample error
TrainErr.elnet = ComputeTrainingError(d.training, BestCoefficients.elnet)
reconstruction = ReconstructionOfTrainingSet(d.training, BestCoefficients.elnet)
###### Now check the quality of the inference
rho.elnet = jacobiano_analitico(JacobianName,
                                BestCoefficients.elnet$J, ncol(d), t.min, t.min+length.training-1)
elnet.matrix.inference = rho.elnet$correlation.matrix
##### Only here we can take the test set out so that model selection did not use the test data
training.error.rho = TrainErr.elnet$rho
training.error.rmse = compute.rmse.train(d.training,reconstruction)
test.error.rho = MeanCorrelation(d.testing, prd)$correlation
test.error.rmse = MeanCorrelation(d.testing, prd)$rmse

##### Make a couple of plots to visualize the results
if(ShowPlot == TRUE){
  source('PlotFunctions.r')
  par(mfrow = c(ncol(d.training), 1))
  for(i in 1:ncol(d.training)){
    Put_train_test_together(d.training[2:nrow(d.training),],d.testing,
                            reconstruction, prd, i)
  }
  par(mfrow = c(1,1))
}
##### Print key information
cat('Out of sample correlation coefficient:', test.error.rho, '\n')
cat('Correlation matrix of the jacobian coefficients:\n')
print(elnet.matrix.inference)
cat('Mean correlation coefficients of the Jacobian:', 
    mean(elnet.matrix.inference[!is.na(elnet.matrix.inference)]),'\n')
cat('Best kernel:', Kernel.Options[other.parameters[idx,1]], '\n')
cat('Best alpha:', other.parameters[idx,1], '\n')
