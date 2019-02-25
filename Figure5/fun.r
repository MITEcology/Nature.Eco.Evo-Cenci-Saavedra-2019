find_optimal_embedding <- function(X,Y, ts, siz.lib, n.samp){
  ### Find the optimal embedding dimension for convergent cross-mapping
  ### Input:
  ### X = dest
  ### Y = pred.spec
  ###
  ### Output:
  ### optimal embedding dimension
  ts.sub = ts[,c(X,Y)]
  ts.sub = scale(ts.sub)
  ###
  loop_cycle <- function(i, X,Y,ts.sub, siz.lib,n.samp){
    dest_xmap_pred <- ccm(as.data.frame(ts.sub), E = i, lib_column = X, 
                          target_column = Y, lib_sizes = siz.lib, num_samples = n.samp, 
                          random_libs = TRUE, replace = TRUE)
    pred_xmap_dest <- ccm(as.data.frame(ts.sub), E = i, lib_column = Y, 
                          target_column = X, lib_sizes = siz.lib, num_samples = n.samp, 
                          random_libs = TRUE, replace = TRUE)
  
    t_xmap_pr_means <- ccm_means(dest_xmap_pred)
    pr_xmap_t_means <- ccm_means(pred_xmap_dest)
  
    y1 <- pmax(0, t_xmap_pr_means$rho)
    y2 <- pmax(0, pr_xmap_t_means$rho)
    rho = max(y1[length(y1)], y2[length(y2)])
    return(rho)
  }
  #### Search in parallel for the best embedding
  embedding = seq(5,25,1)
  Lavoratori = detectCores() - 2
  cl <- makeCluster(Lavoratori, type = "FORK")
  cross.mapping.skills = unlist(parLapply(cl, embedding, loop_cycle, X, Y,ts.sub, siz.lib,n.samp))
  opt.emb = which(cross.mapping.skills == max(cross.mapping.skills))
  stopCluster(cl)
  return(list( optimum.embeddig.dimension = opt.emb, embedding.path = cross.mapping.skills))
}

create_time_series_population <- function(file.to.read, surr = FALSE, which.surrogate, ...){
  d = as.matrix(read.table(file.to.read, header = T))
  if(surr == TRUE){
    number.of.observations = nrow(d)
    number.of.species = ncol(d)
    species.name = colnames(d)
    new_surr = c()
    for (s in 1:ncol(d)){
      if(which.surrogate == 1){
        ### Seasonal surrogate data
        s1 = make_surrogate_data(d[,s], method = "seasonal", num_surr = 1, T_period = 12)
      }else if (which.surrogate == 2){
        ### aaft
        s1 = surrogate(d[,s], method = 2)
      }
      new_surr = cbind(new_surr,s1)
    }
    d = matrix(as.numeric(new_surr), number.of.observations, number.of.species)
    colnames(d) = species.name
  }
  return(d)
}


surrogate_time_series_causality <- function(file.to.read, number_surr, which.surrogate, emb, n.samp,
                                            dest, pred.spec,siz.lib){

  y1_ = y2_ = list()
  for(surr_ in 1:number_surr){
    surr_d = create_time_series_population(file.to.read, surr = TRUE, which.surrogate)
    ts.sub = surr_d[,c(dest,pred.spec)]
    ts.sub = scale(ts.sub)

    dest_xmap_pred <- ccm(as.data.frame(ts.sub), E = emb, lib_column = dest, 
                          target_column = pred.spec, lib_sizes = siz.lib, num_samples = n.samp, 
                          random_libs = TRUE, replace = TRUE)
    pred_xmap_dest <- ccm(as.data.frame(ts.sub), E = emb, lib_column = pred.spec, 
                          target_column = dest, lib_sizes = siz.lib, num_samples = n.samp, 
                          random_libs = TRUE, replace = TRUE)
    t_xmap_pr_means <- ccm_means(dest_xmap_pred)
    pr_xmap_t_means <- ccm_means(pred_xmap_dest)
    y1_[[surr_]] <- pmax(0, t_xmap_pr_means$rho)
    y2_[[surr_]] <- pmax(0, pr_xmap_t_means$rho)
  }
  return(list(y1_ = y1_, y2_ = y2_))
}
take.quantile <- function(X){
  ### X = surrogate ensemble
  surrogate.vector = t(sapply(X, FUN = cbind))
  new_surrogate.vector = sapply(as.data.frame(surrogate.vector), FUN = quantile, probs = seq(0,1,0.05), na.rm = T)
  return(list(surrogate_05 = new_surrogate.vector[2,], surrogate_95 = new_surrogate.vector[20,]))
}


